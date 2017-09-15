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
The base functions for a pipeline method manager.
"""

import os
import atexit
import shutil
import tempfile
from datetime import datetime

try:
    from django.conf import settings
except ImportError:
    settings = {}

try:
    from django.utils.timezone import make_aware
    from django.utils.timezone import now
except ImportError:
    import pytz
    make_aware = lambda dt: timezone.localize(value, is_dst=None)
    from datetime import datetime
    now = datetime.now


class ManagerBase(object):
    def __init__(self, pipedir=None):
        if pipedir is None:
            self.pipedir = tempfile.mkdtemp(prefix='pipeline-')
            atexit.register(self.clean_up)
        else:
            self.pipedir = pipedir

    def job_fn(self, job_id, ext='pid'):
        """Return the filename of the given job_id and type"""
        return os.path.join(self.pipedir, job_id + '.' + ext)

    def clean_up(self):
        """Deletes all data in the piepline directory."""
        if os.path.isdir(self.pipedir):
            shutil.rmtree(self.pipedir)

    def job_read(self, job_id, ext='pid'):
        """Returns the content of the specific job file"""
        fn = self.job_fn(job_id, ext)
        if os.path.isfile(fn):
            with open(fn, 'r') as fhl:
                dt = datetime.fromtimestamp(os.path.getmtime(fn))
                return (make_aware(dt), fhl.read().strip())
        else:
            return (None, None)

    def job_clean(self, job_id, ext):
        """Delete files once finished with them"""
        fn = self.job_fn(job_id, ext)
        if os.path.isfile(fn):
            os.unlink(fn)

    def job_write(self, job_id, ext, data):
        """Write the data to the given job_id record"""
        fn = self.job_fn(job_id, ext)
        with open(fn, 'w') as fhl:
            fhl.write(str(data))


    @property
    def name(self):
        return type(self).__module__.split('.')[-1]

