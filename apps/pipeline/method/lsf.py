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
This is the LSF process for submitting jobs to a compute cluster.

This is used by orchestra's compute cloud.
"""

import os
import shutil
import logging

from datetime import datetime
from subprocess import Popen, PIPE
from django.conf import settings

from .base import ManagerBase

def which(file):
    """In python3.3+ this can be replaced"""
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, file)):
            return os.path.join(path, file)


class JobManager(ManagerBase):
    def __init__(self, *args, **kw):
        for prog in ('bsub', 'bjobs'):
            if not which(prog):
                logging.warn("%s program is not available!" % prog)

        self.group = getattr(settings, 'PIPELINE_LSF_GROUP', 'pipeline')
        self.queue = getattr(settings, 'PIPELINE_LSF_QUEUE', 'short')

        super(JobManager, self).__init__(*args, **kw)

    def submit(self, job_id, cmd, depends=None):
        """
        Open the command locally using bash shell.
        """
        extra = []
        if depends:
            extra += ['-w', 'done(%s)' % depends]
        Popen(['bsub', '-J', job_id, '-g', self.group, '-q', self.queue] + \
              extra + ['-o', self.job_fn(job_id, 'err'), '-W', '2:00', cmd],
            shell=False, stdout=None, stderr=None, close_fds=True)
        return True

    def stop(self, job_id):
        """Stop the given process using bkill"""
        Popen(['bkill', '-J', job_id, '-g', self.group])

    def status(self, job_id, clean=False):
        """Returns if the job is running, how long it took or is taking and other details."""
        # Get the status for the listed job, how long it took and everything
        p = Popen(['bjobs', '-J', job_id, '-a', '-W'], stdout=PIPE, stderr=None)
        (out, err) = p.communicate()

        #JOBID      USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME  PROJ_NAME CPU_USED MEM SWAP PIDS START_TIME FINISH_TIME
        #558531  mo131   RUN   short      gallery.orchestra ottavino002-229.orchestra sleep 50   02/03-14:55:11 default    000:00:00.00 2160   230600 5866,5869,5873 02/03-14:55:19 - 

        # Prevent spaces in the job_id from hurting the processing
        out = out.replace(job_id, 'job_id_here')

        # Turn the output into a dictionary useful
        lines = out.split('\n')
        if len(lines) <= 1:
            return {}

        data = zip(lines[0].lower().split(), lines[1].split())

        # When getting date-times from lsf, we have to convert and add the year (very odd)
        year = "%d/" % datetime.now().year
        for dkey in data:
            if not dkey.endswith('_time'):
                continue
            if data[dkey] == '-':
                data[dkey] = None
            else:
                data[dkey] = datetime.strptime(year + data[dkey], '%Y/%m/%d-%H:%M:%S')

        status = {
            'RUN': 'running', # Active
            'PEND': 'pending', # Stopped because we asked it to be
            'PSUSP': 'finished',
            'SSUSP': 'finished',
            'DONE': 'finished', # Pining for the fyords
            'EXIT': 'finished',
        }.get(data['stat'])

        ret = 0 if data['stat'] == 'DONE' else 1
        (_, err) = self.job_read(job_id, 'err')

        if status == 'finished' and clean:
            self.job_clean(job_id, 'err')

        return {
            'submitted': data['submit_time'],
            'started': data['start_time'],
            'finished': data['finish_time'],
            'pid': data['JOBID'],
            'status': status,
            'return': ret,
            'error': err,
        }
        

