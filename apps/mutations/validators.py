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
Provide low level validators for database fields. Provides octal validation.
"""

from django.core.exceptions import ValidationError
from django.utils.translation import ugettext_lazy as _

def is_octal(value):
    """Validate that an incoming value is an octal value"""
    value = str(value)
    if not value.isdigit() or "8" in value or "9" in value:
        raise ValidationError(_('%(value)s is not an octal number'), params={'value': value})

