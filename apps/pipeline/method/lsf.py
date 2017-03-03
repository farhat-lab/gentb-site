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

from datetime import datetime

from subprocess import Popen, PIPE

from .base import ManagerBase

LSF_QUEUE = 'short'

class Manager(ManagerBase):
    def __init__(self):
        # XXX check is bsub and bjobs are available commands
        # TODO check is bsub has queues and if it has the configured queu
        self.group = 'pipeline' # TODO Get this from the django settings.
        self.ready = True

    def submit(self, job_id, cmd):
        """
        Open the command locally using bash shell.
        """
        Popen(['bsub', '-J', job_id, '-g', self.group, '-q', QUEUE,
            '-W', '2:00', cmd],
            shell=False, stdout=None, stderr=None, close_fds=True)

    def stop(self, job_id):
        """Stop the given process using bkill"""
        Popen(['bkill', '-J', job_id, '-g', self.group])

    def status(self, job_id):
        """Returns if the job is running, how long it took or is taking and other details."""
        """
        started - datetime
        finished - datetime
        pid - integer
        job_id - string
        status - string / fixed list
        return code - number
        error - multi-line string (stderr)
        """

        # Get the status for the listed job, how long it took and everything
        p = Popen(['bjobs', '-J', job_id, '-a', '-W'], stdout=PIPE, stderr=None)
        (out, err) = p.communicate()

        #JOBID      USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME  PROJ_NAME CPU_USED MEM SWAP PIDS START_TIME FINISH_TIME
        #558531  mo131   RUN   short      gallery.orchestra ottavino002-229.orchestra sleep 50   02/03-14:55:11 default    000:00:00.00 2160   230600 5866,5869,5873 02/03-14:55:19 - 

        # Prevent spaces in the job_id from hurting the processing
        out = out.replace(job_id, 'job_id_here')

        # Turn the output into a dictionary useful
        lines = out.split('\n')
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

        data['running'] = data.pop('stat') == 'RUN'

        return data
        

