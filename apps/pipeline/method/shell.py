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
This is the most basic way of running jobs, in the local shell.
"""

import os
import shutil
import signal

from subprocess import Popen
from datetime import datetime

from .base import ManagerBase

class JobManager(ManagerBase):
    """
    The job is submitted to the raw system and is controlled via it's pid.
    This pid must be stored in the pipeline-shell directory in tmp.
    """
    def clean_up(self):
        """Deletes all data in the piepline directory."""
        if os.path.isdir(self.pipedir):
            shutil.rmtree(self.pipedir)

    def job_read(self, job_id, ext='pid'):
        """Returns the content of the specific job file"""
        fn = self.job_fn(job_id, ext)
        if os.path.isfile(fn):
            with open(fn, 'r') as fhl:
                return (datetime.fromtimestamp(os.path.getmtime(fn)),
                        fhl.read().strip())
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

    watch = os.path.join(os.path.dirname(__file__), 'watch.py')
    DEP = watch + ' %(fn)s && '
    RET = '; echo $? > "%(fn)s"'
    def submit(self, job_id, cmd, depends=None):
        """
        Open the command locally using bash shell.
        """
        if depends:
            (_, pid) = self.job_read(depends, 'pid')
            (_, ret) = self.job_read(depends, 'ret')
            if pid:
                if ret not in (0, None):
                    # Refuse to submit job if the dependant job failed
                    return False
                elif ret is None:
                    # If we depend on another process, then watch for it's
                    # return by watching for the appearence of the ret file
                    # then checking it's content
                    cmd = (self.DEP % {'fn': self.job_fn(depends, 'ret')}) + cmd

        # Collect the standard error into an err file
        err = open(self.job_fn(job_id, 'err'), 'w')
        err.write('S')
        err.seek(0)
        # Dump the standard out into oblivion
        out = open('/dev/null', 'w')
        # Collect the return code into the ret file
        cmd += self.RET % {'fn': self.job_fn(job_id, 'ret')}

        # Run the large shell command
        proc = Popen(cmd, shell=True, stdout=out, stderr=err, close_fds=True)

        # Pid the pid (process id) into the pid file
        self.job_write(job_id, 'pid', proc.pid)
        return True

    def all_children(self):
        """Yields all running children remaining"""
        from collections import defaultdict
        pids = defaultdict(list)
        for pid in os.listdir('/proc'):
            if pid.isdigit() and self.is_running(pid):
                pid_fn = "/proc/%s/status" % str(pid)
                data = dict(line.split(':\t', 1) for line in open(pid_fn).readlines())
                pids[int(data['PPid'])].append(int(pid))
        parents = [os.getpid()]
        while parents:
            for child in pids.get(parents.pop(0), []):
                if self.state_and_clear(pid, False):
                    parents.append(child)
                    yield child
                    
    def clean_up(self):
        """Create a list of all processes and kills them all"""
        pids = list(self.all_children())
        for pid in pids:
            try:
                os.kill(int(pid), signal.SIGKILL)
                os.waitpid(int(pid), 0)
            except OSError:
                pass
        return pids

    def stop(self, job_id):
        """Send a SIGTERM to the job and clean up"""
        (_, pid) = self.job_read(job_id, 'pid')
        if pid is None:
            return
        if self.is_running(pid):
            try:
                os.kill(int(pid), signal.SIGKILL)
                (_, ret) = os.waitpid(int(pid), 0)
                self.job_write(job_id, 'err', 'S')
                self.job_write(job_id, 'ret', ret)
            except IOError:
                pass

    def is_running(self, pid):
        """Returns true if the process is still running"""
        return os.path.exists("/proc/%d/status" % int(pid))

    def status(self, job_id, clean=False):
        """Returns a dictionary containing status information,
        can only be called once as it will clean up status files!"""
        (started, pid) = self.job_read(job_id, 'pid')
        (finished, ret) = self.job_read(job_id, 'ret')
        (_, err) = self.job_read(job_id, 'err')
        status = 'finished'
        if err == 'S':
            err = ''
            if ret != '0':
                status = 'stopped'

        if pid is None and ret is None:
            return None

        status = self.state_and_clear(pid, status)

        if clean and (finished or err):
            self.job_clean(job_id, 'pid')
            self.job_clean(job_id, 'ret')
            self.job_clean(job_id, 'err')

        return {
          'status': status,
          'started': started,
          'finished': finished,
          'return': int(ret) if ret is not None else None,
          'error': err,
        }

    def state_and_clear(self, pid, default=None):
        """Gets the status of the process and waits for zombies to clear"""
        pid_fn = "/proc/%d/status" % int(pid)
        if os.path.exists(pid_fn):
            data = dict(line.split(':\t', 1) for line in open(pid_fn).readlines())
            if data['State'][0] == 'Z':
                # Clear up zombie processes waiting for us to pid them.
                os.waitpid(int(pid), 0)
                return default
            return {
              'D': 'sleeping', # Machine is too busy
              'S': 'sleeping', # Busy doing nothing
              'R': 'running', # Active
              'T': 'pending', # Stopped because we asked it to be
              'X': 'finished', # Pining for the fyords
            }[data['State'][0]]
        return default

