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
This allows tests to operate pipelines without doing anything.
"""
from django.utils.timezone import now
from .base import ManagerBase

class JobManager(ManagerBase):
    cmds = {}

    def submit(self, job_id, cmd, depends=None):
        if 'bad-submit' in job_id:
            return False
        # Command, sleeping, started, completed, err
        self.cmds[job_id] = (cmd, False, False, False, False)
        return True

    def run_all(self):
        for job_id in self.cmds:
            self.cmds[job_id][2] = True

    def finish_all(self):
        for job_id in self.cmds:
            self.cmds[job_id][3] = True

    def clean_up(self):
        self.cmds = {}

    def stop(self, job_id):
        if self.cmds[job_id][2] and not self.cmds[job_id][3]:
            self.cmds[job_id][3] = True
            self.cmds[job_id][4] = True

    def status(self, job_id, clean=False):
        if job_id in self.cmds:
            cmd = self.cmds[job_id]
        else:
            def parse(v):
                try:
                    return int(v)
                except:
                    return v
            cmd = [parse(i) for i in job_id.split('_')]
            cmd = cmd + [False, False, False, False]

        status = 'pending'
        if cmd[4]:
            status = 'finished' # Error
        elif cmd[3]:
            status = 'finished' # Non-error
        elif cmd[1]:
            status = 'sleeping'
        elif cmd[2]:
            status = 'running'

        return {
          'status': 'finished',
          'started': now() if cmd[2] else None,
          'finished': now() if cmd[3] else None,
          'return': int(cmd[4]),
          'error': '',
        }

