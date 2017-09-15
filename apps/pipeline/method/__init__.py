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

import inspect

from importlib import import_module

from .base import ManagerBase, settings

DEFAULT = 'apps.pipeline.method.shell'

def get_job_manager(module_id=None, pipe_root=None):
    if pipe_root is None:
        pipe_root = getattr(settings, 'PIPELINE_ROOT', None)
    if module_id is None:
        module_id = getattr(settings, 'PIPELINE_MODULE', DEFAULT)

    # Already a job manager, so return
    if isinstance(module_id, ManagerBase):
        return module_id

    # A job manager class, create object and return
    if inspect.isclass(module_id):
        return module_id(pipedir=pipe_root)

    # A name to a job manage, import and create
    try:
        module = import_module(module_id)
        obj = getattr(module, 'JobManager')(pipedir=pipe_root)
        return  obj
    except Exception as err:
        raise ValueError("Pipeline module is not a pipeline method or "
                "not configured correctly, %s: %s" % (module_id, str(err)))

