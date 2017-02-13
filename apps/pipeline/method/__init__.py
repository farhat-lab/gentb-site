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
Gets the configured pipeline module and initialises it.
"""

from importlib import import_module
from django.conf import settings

module_id = getattr(settings, 'PIPELINE_MODULE', 'apps.pipeline.method.shell')

try:
    module = import_module(module_id)
    JobManager = getattr(module, 'JobManager')()
except Exception as err:
    raise ValueError("Pipeline module is not a pipeline method or "
            "not configured correctly, %s: %s" % (module_id, str(err)))

