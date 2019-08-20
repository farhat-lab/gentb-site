#
# Copyright (C) 2018 Maha Farhat
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
Common useful pipeline utilities.
"""

from collections import defaultdict

def file_as_inputs(m2m, save_to=None):
    """Adds each m2m to the save_to dictionary"""
    if not save_to:
        save_to = defaultdict(list)
    for pfa in m2m.all():
        try:
            save_to[pfa.name].append(str(pfa.store.file))
        except IOError:
            save_to[pfa.name].append('XX:%s' % pfa.store.name)
    return save_to

