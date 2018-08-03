#
# Copyright (C) 2015 Adam Bogdal
#               2017 Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Code from: https://github.com/bogdal/django-dropboxchooser-field by Adam Bogdal.
# It is unlicensed. We're assuming Public Domain, MIT or some other permissive
# license.
#

from django.core.urlresolvers import reverse
from django.conf import settings
from django.forms import TextInput, Textarea

import json

class UploadChooserWidget(TextInput):
    input_type = 'upload-chooser'

    def __init__(self, extensions=None, attrs=None, buckets=None):
        extensions = ['.' + x.strip().strip('.') for x in (extensions or [])]
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
        attrs.pop('required', None)
        attrs['data-resumable_url'] = reverse('uploads:resumable')
        attrs['data-manual_url'] = reverse('uploads:manual')
        ret = render(name, value, attrs)
        for bucket, match, label, link in self.buckets:
            # This is evil, but django changed how it works
            self.input_type = 'bucket'
            kw = attrs.copy()
            kw['data-match'] = match
            kw['data-parent'] = kw['id']
            kw['data-label'] = label
            kw['data-bucket'] = bucket
            kw['data-link'] = link
            kw['id'] += '_' + bucket
            ret += render(name + '_' + bucket, '', kw)
            self.input_type = 'upload-chooser'
        return ret

    class Media:
        js = [
          # Dropbox javascript support
          'https://www.dropbox.com/static/api/2/dropins.js?cache=2',
          # Chunked file uploader support
          'js/resumable.js',
          # Bootbox used in URL Uploader
          'js/bootbox.min.js',
          # Generic uploader support (binds together the above)
          'js/uploads/chooser.js']
        css = {'all': ('css/uploads/chooser.css',)}


class UploadTableWidget(Textarea):
    """
    Show an importer UI which allows a file to be selected from the user's
    hard disk and then that file is parsed on the client side to check
    that it matches the definition or other checks.
    """
    def __init__(self, columns, parsers=None):
        kw = { 
            'style': 'display: none',
            'data-columns': ','.join(columns),
            'type': 'upload-table',
        }

        # Convert list to dictionary using the columns as a guide.
        if isinstance(parsers, (list, tuple)):
            parsers = dict(zip(columns[:len(parsers)], parsers))

        # Add each parser as it's own data attribute on the Text element.
        for key, parser in (parsers or {}).items():
            if parser is not None:
		kw['data-parser_' + key] = unicode(parser)

        super(UploadTableWidget, self).__init__(kw)


    class Media:
        js = [
            'js/papaparse.min.js',
            'js/uploads/importer.js',
        ]
        css = {'all': ('css/uploads/importer.css',)}

