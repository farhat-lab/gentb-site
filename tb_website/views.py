#
# Copyright 2016, Martin Owens <doctormo@gmail.com>
#           2017, Maha Farhat for GenTB
#
# This file is part of the software inkscape-web, consisting of custom 
# code for the Inkscape project's django-based website.
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
"""
Template Error
"""

from django.views.generic import TemplateView

class Error(TemplateView):
    @classmethod
    def as_error(cls, status):
        view = cls.as_view(template_name='error/%s.html' % status)
        def _inner(request):
            response = view(request, status=int(status))
            if hasattr(response, 'render'):
                response.render()
            return response
        return _inner

    def post(self, request, **kw):
        return self.get(request, **kw)

    def get(self, request, **kw):
        context = self.get_context_data(**kw)
        return self.render_to_response(context, **kw)


