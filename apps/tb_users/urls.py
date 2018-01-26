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
Login and user views which use the 'users:' prefix.
"""

from django.conf.urls import url
from django.contrib.auth.views import LoginView, LogoutView

from .views import UserSignUp, SignUpSuccess

urlpatterns = [
  url(r'^signup/$', UserSignUp.as_view(), name="signup"),
  url(r'^signup/success/$', SignUpSuccess.as_view(), name="signup_success"),

  url(r'^login/', LoginView.as_view(), name='login'),
  url(r'^logout/', LogoutView.as_view(), name='logout'),
]

