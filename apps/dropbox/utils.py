
import logging
import requests

from os.path import join, getsize
from django.core.files.storage import FileSystemStorage

class Download(object):
    """Wrap the requests module for django's storage backend."""
    def __init__(self, url):
        self.filepath = None
        self.io = requests.get(url, stream=True)

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

