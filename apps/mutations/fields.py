#
# Copyright (C) 2016   Dr. Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the 
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Input fields for manual mutation entry.
"""

from django.forms.fields import CharField
from django.forms.widgets import Textarea

class ListerWidget(Textarea):
    class Media:
        js = ('js/lister.js',)

class GeneticInputField(CharField):
    def __init__(self, data_url, *args, **kw):
        self.data_url = data_url
        kw['widget'] = ListerWidget
        super(GeneticInputField, self).__init__(*args, **kw)

    def to_python(self, value):
        "Returns a Unicode object."
        if value in self.empty_values:
            return ''
        return value

    def widget_attrs(self, widget):
        attrs = super(CharField, self).widget_attrs(widget)
        attrs['data-data-url'] = self.data_url
        attrs['class'] = 'lister'
        return attrs

