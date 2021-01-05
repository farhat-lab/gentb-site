#
# Copyright 2015, Martin Owens
#           2017, MIT TESS
#           2018, Maha Farhat
#
# tess-tev is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tess-tev is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with tess-tev.  If not, see <http://www.gnu.org/licenses/>.
#
# Code imported from inkscape-web 2017-10-18 AGPLv3
# Code imported from tev.mit.edu 2018-09-31 AGPLv3
#

from django.contrib.auth.views import redirect_to_login
from django.core.exceptions import PermissionDenied

class ProtectedMixin():
    """Combine login required and permission specific as one mixin"""
    def is_permitted(self, user):
        """Is the user allowed to access this resource"""
        if hasattr(self, 'permission'):
            if user is None or not user.has_perm(self.permission):
                return self.not_allowed()
        elif hasattr(self, 'permissions'):
            if user is None or not user.has_perms(self.permissions):
                return self.not_allowed()
        elif getattr(self, 'super_only', False):
            if user is None or not user.is_superuser:
                return self.not_allowed()
        elif getattr(self, 'staff_only', False):
            if user is None or not user.is_staff:
                return self.not_allowed()
        return True

    @staticmethod
    def not_allowed():
        """Action to take if the user has the wrong permissions"""
        raise PermissionDenied

    def not_authenticated(self):
        """Action to take if the user isn't identified"""
        return redirect_to_login(self.request.build_absolute_uri())

    def dispatch(self, request, *args, **kw):
        """When this view is used, protect with authorisation steps"""
        user = request.user if request.user.is_authenticated else None
        if not self.is_permitted(user):
            if user is None:
                return self.not_authenticated()
            return self.not_allowed()
        return super(ProtectedMixin, self).dispatch(request, *args, **kw)
