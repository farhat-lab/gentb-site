from os.path import isdir, join
import urllib2
import shutil
from datetime import datetime

class DropBoxFileDownloader:
    """
    Given a Dropbox Public link:
    - Download the related files in .zip format
    - Unzip
    - Keep only the VCF and FastQ files
    """
    def __init__(self, dropbox_link, target_directory):
        if dropbox_link is None:
            raise ValueError("dropbox_link cannot be None")
        if not isdir(target_directory):
            raise IOError("target_directory not found: %s" % target_directory)

        self.dropbox_link = dropbox_link
        self.target_directory = target_directory
        self.err_found = False
        self.err_msg = None

        self.download_files()

    def add_err(self, m):
        self.err_found = True
        self.err_msg = m

    def is_valid_dropbox_link(self):
        d = self.dropbox_link.lower()

        dbox_start = 'https://www.dropbox.com'
        if not d.startswith(dbox_start):
            self.add_err("Not a dropbox url.  Doesn't start with '%s'" % dbox_start)
            return False

        if d.endswith('?dl=0'):
            self.dropbox_link = self.dropbox_link.replace('?dl=0', '?dl=1')
            return True

        if d.endswith('?dl=1'):
            return True

        self.add_err('Does not appear to be a valid download link.  Should end with "?dl=0" or "?dl=1"')
        return False

    def download_files(self):
        print 'download_files'
        if not self.is_valid_dropbox_link():
            return

        dbox_zip_fname = 'zippity.zip'
        dbox_zip_fullname = join(self.target_directory, dbox_zip_fname)
        chunk_size = 16 * 1024
        req = urllib2.urlopen(self.dropbox_link)
        print 'download link: {0}'.format(self.dropbox_link)

        print 'pre download: {0}'.format(datetime.now())
        with open(dbox_zip_fullname, 'wb') as fp:
            #cnt += 1
            print 'downloading...'
            shutil.copyfileobj(req, fp, chunk_size)

        print 'file written: {0}'.format(dbox_zip_fullname)
        print 'post download: {0}'.format(datetime.now())

if __name__ == '__main__':
    # 1k
    #dlink = 'https://www.dropbox.com/s/xxlpguql6gtq311/001.txt?dl=0'
    # 2mb
    #dlink = 'https://www.dropbox.com/s/omnkw09mdnd1r51/2mb.txt?dl=0'
    # 1gb
    #dlink = 'https://www.dropbox.com/s/ls4xqbyvd3lj3p9/1gb.txt?dl=0'
    dlink = 'https://www.dropbox.com/sh/2rr42u7v7x9549x/AAAcg8v8AzmWva0Nx-kwr7qxa?dl=0'
    dbox = DropBoxFileDownloader(dlink, 'downloaded')

"""
from dropbox_info import *

link = 'https://www.dropbox.com/sh/r7qb9skx1vc3xxd/AABPxsafyiXtcsJ0Vg20JNloa?dl=0'
print """curl -X POST https://api.dropbox.com/1/metadata/link -u {0}:{1} -d link={2}""".format(app_key, secret_key, link)

app curl -X POST https://api.dropbox.com/1/metadata/link -u <APP_KEY>:<APP_SECRET> \
      -d link=https://www.dropbox.com/sh/748f94925f0gesq/AAAMSoRJyhJFfkupnAU0wXuva?dl=0
"""
