#
# Copyright (C) 2017 Maha Farhat
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
"""
The actual input, formatting and etc.
"""
import json
from collections import defaultdict

from django.forms.fields import Field

from .widgets import UploadChooserWidget
from .models import UPLOADERS

class UploadField(Field):
    widget = UploadChooserWidget

    def __init__(self, **kw):
        self.directory = kw.pop('dir', None)
        self.extensions = kw.get('extensions', None)
        if 'widget' not in kw:
	    kw['widget'] = UploadChooserWidget(
                extensions=kw.pop('extensions', None),
                attrs=kw.pop('attrs', None),
                buckets=kw.pop('buckets', None))
        super(UploadField, self).__init__(**kw)

    def to_python(self, value):
        """Returns the correct upload file type."""
        try:
            raw = json.loads(value)
        except:
            raw = []

        ret = defaultdict(list)
        for datum in raw:
            prefix = self.get_prefix(datum['id'])
            cls = UPLOADERS[datum.get('source', '')]
            bucket = datum['bucket']
            ret[bucket].append(cls.build_upload(prefix, datum))
        return ret

    def get_prefix(self, filename):
        for ext in self.extensions or []:
            if filename.endswith(ext):
                return filename[:len(filename)-len(ext)]
        return filename

