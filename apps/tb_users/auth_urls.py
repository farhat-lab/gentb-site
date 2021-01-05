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

from django.conf.urls import url, include
from django.contrib.auth.views import (
  PasswordResetView, PasswordResetConfirmView,
  PasswordResetCompleteView, PasswordResetDoneView,
)

def url_tree(regex, *urls):
    class UrlTwig():
        urlpatterns = urls
    return url(regex, include(UrlTwig))

UIDB = r'^(?P<uidb64>[0-9A-Za-z_\-]+?)/(?P<token>.+)/$'

urlpatterns = [
  url_tree(r'^password_reset/',
    url(r'^$',      PasswordResetView.as_view(), name='password_reset'),
    url(UIDB,       PasswordResetConfirmView.as_view(), name='password_reset_confirm'), 
    url(r'^done/$', PasswordResetCompleteView.as_view(), name='password_reset_complete'),
    url(r'^sent/$', PasswordResetDoneView.as_view(), name='password_reset_done'),
  )
]

