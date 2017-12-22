#
# Copyright 2015, Martin Owens
#
# inkscape-web is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# inkscape-web is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with inkscape-web.  If not, see <http://www.gnu.org/licenses/>.
#

from datetime import datetime

try:
    from django.template.library import Library
except ImportError:
    # Django <= 1.9
    from django.template import Library

register = Library()

@register.filter("epoch")
def epoch(seconds):
    return datetime.fromtimestamp(seconds)

