
import uuid
import logging
import ftplib
import requests
import socket

from os.path import join, getsize
from django.core.files.storage import FileSystemStorage

def get_uuid():
    return uuid.uuid4().hex

class Download(object):
    """Wrap the requests module for django's storage backend."""
    def __init__(self, url):
        self.filepath = None
        self.io = requests.get(url, stream=True)

    def is_ok(self):
        return self.io.status_code == 200

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
    # XXX We could store the username/password and strip it out
    # of the returned value to make it easier for us to control
    # the ftp password (and not store it a lot of times)
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
    except Exception: # XXX Replace with specifics.
        pass # Get parent of path and try again.

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


