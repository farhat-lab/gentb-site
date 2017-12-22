#
# Copyright 2017, Maha Farhat
#
# versioner is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# versioner is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with versioner.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Show any incoming versions to the website user.
"""

from django.views.generic import TemplateView
from django.utils.decorators import method_decorator
from django.contrib.auth.decorators import login_required

class Version(TemplateView):
    template_name = 'versioner/report.html'

    @method_decorator(login_required)
    def dispatch(self, request, *args, **kwargs):
        """Only allow a logged in users to view"""
        return super(Version, self).dispatch(request, *args, **kwargs)

