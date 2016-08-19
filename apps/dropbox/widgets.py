#
# This code original a part of https://github.com/bogdal/django-dropboxchooser-field
# which is unlicensed. We're assuming Public Domain, MIT or some other permissive
# license.
#
# Heavily modified.
#

from django.conf import settings
from django.forms import TextInput

import json

class DropboxChooserWidget(TextInput):
    input_type = 'dropbox-chooser'

    def __init__(self, extensions=None, attrs=None):
        extensions = ['.' + x.strip('.') for x in (extensions or [])]
        kw = {
            'style': 'display: none',
            'data-app-key': getattr(settings, 'DROPBOX_APP_KEY', None),
            'data-extensions': " ".join(extensions),
        }
        kw.update(attrs or {})
        super(DropboxChooserWidget, self).__init__(kw)

    class Media:
        js = ['https://www.dropbox.com/static/api/2/dropins.js?cache=2',
              'js/dropbox/chooser.js?cache=2']
        css = {'all': ('css/dropbox/chooser.css?cache=2',)}

