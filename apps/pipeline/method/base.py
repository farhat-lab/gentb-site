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
import sys
import atexit
import shutil
import tempfile
from datetime import datetime

try:
    from django.conf import settings #pylint: disable=unused-import
except ImportError:
    settings = {}

try:
    from django.utils.timezone import make_aware
    from django.utils.timezone import now #pylint: disable=unused-import
except ImportError:
    import timezone
    make_aware = lambda dt: timezone.localize(dt, is_dst=None)
    now = datetime.now


class ManagerBase(object):
    """Manage any number of pipeline methods such as shell, slurm, lsb, etc"""
    name = property(lambda self: type(self).__module__.split('.')[-1])

    def __init__(self, pipedir=None):
        if pipedir is None:
            self.pipedir = tempfile.mkdtemp(prefix='pipeline-')
            atexit.register(self.clean_up)
        else:
            self.pipedir = pipedir

    def job_fn(self, job_id, ext='pid'):
        """Return the filename of the given job_id and type"""
        if not os.path.isdir(self.pipedir):
            os.makedirs(self.pipedir)
        return os.path.join(self.pipedir, job_id + '.' + ext)

    def clean_up(self):
        """Deletes all data in the piepline directory."""
        if os.path.isdir(self.pipedir):
            shutil.rmtree(self.pipedir)

    def job_read(self, job_id, ext='pid'):
        """Returns the content of the specific job file"""
        filen = self.job_fn(job_id, ext)
        if os.path.isfile(filen):
            with open(filen, 'r') as fhl:
                dtm = datetime.fromtimestamp(os.path.getmtime(filen))
                return (make_aware(dtm), fhl.read().strip())
        else:
            return (None, None)

    def job_clean(self, job_id, ext):
        """Delete files once finished with them"""
        filen = self.job_fn(job_id, ext)
        if os.path.isfile(filen):
            os.unlink(filen)
            return True
        return False

    def job_write(self, job_id, ext, data):
        """Write the data to the given job_id record"""
        filen = self.job_fn(job_id, ext)
        with open(filen, 'w') as fhl:
            fhl.write(str(data))

    def job_stale(self, job_id):
        """Figure out if a job has stale return files"""
        if self.job_clean(job_id, 'ret'):
            sys.stderr.write("Stale job file cleared: {}\n".format(job_id))
            self.job_clean(job_id, 'pid')
