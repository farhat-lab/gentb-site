
import json
import uuid
import socket
import ftplib
import logging
import requests
import inspect

from md5 import md5
from os.path import join, getsize
from django.conf import settings
from django.core.files.storage import FileSystemStorage

class CachedUrl(object):
    """Saves a URL and replaces it with a hash that can be loaded server side"""
    def __init__(self, url):
        if url.startswith('url://'):
            self.name = url.split('://')[-1].split('/')[0]
        else:
            self.name = md5(url).hexdigest()
            (self.protocol, self.url) = url.split('://', 1)

        self.loc = os.path.join(settings.UPLOAD_CACHE_ROOT, self.name + '.json')

        if os.path.isfile(self.loc):
            (url_hash, _) = self.url.split('://')[-1].split('/', 1)
            with open(self.url_cache_filename(url_hash)) as fhl:
                data = json.loads(fhl.read())
                self.protocol = data['protocol']
                self.url = data['url']

    def save(self):
        """
        Saves the url and protocol into a server side cache, so
        details about a file's location isn't beamed around all the time.
        """
        with open(self.loc, 'w') as fhl:
            fhl.write(json.dumps(dict(protocol=self.protocol, url=self.url)))

    def url(self, fn):
        if self.url.endswith(fn):
            self.url = self.url[:-len(fn)]
        self.save()
        return "url://%s" % os.path.join(self.url, fn)


class LocalFile(object):
    """Fake download for a local file already on the server"""
    protocols = ['file']

    def __init__(self, url):
        self.url = url

    def __iter__(self):
        """ 
        Try and load a local directory name, if available.
        """
        # We limit users to only those with this permission
        if not self.request.user.has_perm('uploads.add_manualuploadfile'):
            raise PermissionDenied

        path = self.url
        if os.path.isdir(path):
            for fn in os.listdir(path):
                full = os.path.join(path, fn) 
                if os.path.isfile(full) and self.match_file(fn):
                    yield (fn, os.path.getsize(full))

        elif os.path.isfile(path): # Don't filter with match_file
            yield (os.path.basename(path), os.path.getsize(path))


class DownloadFtp(object):
    protocols = ['ftp', 'ftps']

    def __iter__(self):
        """Load an FTP url"""
        for (name, size) in ftp(url, tls=tls):
            yield (name, size)


def get_uuid():
    return uuid.uuid4().hex


class Download(object): 
    """Wrap the requests module for django's storage backend."""
    protocols = ['http', 'https']

    def __init__(self, url):
        self.filepath = None
        self.io = requests.get(url, stream=True)
        #requests.get(..., auth=(username, password))

    def __iter__(self):
        from bs4 import BeautifulSoup
        # List all the files in a url request.
        with requests.get(self.url, stream=True) as io:
            if io.status_code == 200:
                if io.headers.get('content-type') == 'text/html':
                    page = BeautifulSoup(io.text, 'html.parser')
                    page.prettify()
                    for anchor in page.find_all('a', href=True):
                        yield (anchor['href'], -1)
                else:
                    yield (self.url, io.headers['Content-length'])
            else:
                raise ValueError("Bad http response. ERR:%d" % io.status_code)

    def is_ok(self):
        return self.io.status_code == 200

    def get_error(self):
        return self.io.text

    def save(self, path, filename):
        """Perform the download in chunks"""
        storage = FileSystemStorage(location=path)
        self.filepath = join(path, storage.save(filename, self))

    @property
    def size(self):
        if not self.filepath:
            return 0
        return getsize(self.filepath)

    def chunks(self):
        """Translate django storage backend call, into Response iter"""
        return self.io.iter_content(50 * 1024) # 50KB chunks


def ftp(url, tls=True, **kw):
    """
    Wrapper around ftplib.FTP

     * allows the use of username:password@servername/path full urls
       instead of just server name
     * combines both FTP and FTP_TLS into one class.

    """
    method = [ftplib.FTP, ftplib.FTP_TLS][bool(tls)]
    if '@' in url:
        (kw['user'], url) = url.split('@', 1)
        if ':' in kw['user']:
            (kw['user'], kw['passwd']) = kw['user'].rsplit(':', 1)

    path = '.'
    if '/' in url:
        (url, path) = url.split('/', 1)

    port = 21
    if ':' in url:
        (url, port) = url.rsplit(':', 1)

    if not url:
        raise AttributeError('No ftp url provided')

    try:
        print "CONNECTING TO %s:%s" % (url, str(port))
        ftp = method(timeout=1)
        ftp.connect(url, int(port))
        ftp.login(**kw)
    except socket.gaierror:
        raise AttributeError('Can\'t connect to server "%s"' % url)
    except socket.timeout:
        raise AttributeError('Timeout connecting to server "%s"' % url)

    try:
        ftp.cwd(path)
    except Exception:
        ftp.cmd(os.path.dirname(path))
        filename = os.path.basename(path)

    ret = []
    def _parse_LIST(line):
        """Parses each line from the FTP LIST command"""
        parts = line.split(None, 8)
        if parts[0].startswith('-r'):
            size = int(parts[4])
            name = parts[-1]
            ret.append((name, size))

    ftp.retrlines('LIST', callback=_parse_LIST)
    return ret


DOWNLOADERS = dict()
for cls in locals().values():
    if inspect.isclass(cls) and hasattr(cls, 'protocols'):
        DOWNLOADERS.update([(prot, cls) for prot in cls.protocols])

