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

class UploadChooserWidget(TextInput):
    input_type = 'dropbox-chooser'

    def __init__(self, extensions=None, attrs=None, buckets=None):
        extensions = ['.' + x.strip('.') for x in (extensions or [])]
        kw = {
            'style': 'display: none',
            'data-app-key': getattr(settings, 'DROPBOX_APP_KEY', None),
            'data-extensions': " ".join(extensions),
        }
        self.buckets = buckets or []
        kw.update(attrs or {})
        super(UploadChooserWidget, self).__init__(kw)

    def render(self, name, value, attrs):
        render = super(UploadChooserWidget, self).render
        ret = render(name, value, attrs)
        for bucket, match, label in self.buckets:
            kw = attrs.copy()
            kw['type'] = 'bucket'
            kw['data-extensions'] = match
            kw['data-parent'] = kw['id']
            kw['data-label'] = label
            kw['id'] += '_' + bucket
            ret += render(name + '_' + bucket, '', kw)
        return ret

    class Media:
        js = ['https://www.dropbox.com/static/api/2/dropins.js?cache=2',
              'js/dropbox/chooser.js?cache=2']
        css = {'all': ('css/dropbox/chooser.css?cache=2',)}


# XXX UploadBucket should be it's own hidden widget type.

