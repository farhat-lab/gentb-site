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
from .base import ManagerBase, now

COMMAND, IS_SLEEP, IS_STARTED, IS_COMPLETE, RETURN_CODE, ERROR_OUT = range(6)

class JobManager(ManagerBase):
    cmds = {}

    def submit(self, job_id, cmd, depends=None):
        if 'bad-submit' in job_id:
            return False
        # Command, sleeping, started, completed, err
        self.cmds[job_id] = [cmd, None, None, None, None, '']
        return True

    def run_all(self):
        for job_id in self.cmds:
            self.cmds[job_id][IS_STARTED] = now()

    def finish_all(self, status=0, err=''):
        for job_id in self.cmds:
            self.cmds[job_id][IS_COMPLETE] = now()
            self.cmds[job_id][RETURN_CODE] = status
            self.cmds[job_id][ERROR_OUT] = err

    def clean_up(self):
        self.cmds = {}

    def stop(self, job_id):
        if self.cmds[job_id][IS_STARTED] and not self.cmds[job_id][IS_COMPLETE]:
            self.cmds[job_id][IS_COMPLETE] = True
            self.cmds[job_id][RETURN_CODE] = 127
            self.cmds[job_id][ERROR_OUT] = 'STOPPED'

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
          'status': status,
          'started': cmd[IS_STARTED],
          'finished': cmd[IS_COMPLETE],
          'return': cmd[RETURN_CODE],
          'error': cmd[ERROR_OUT],
        }

