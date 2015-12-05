"""
12/4/2015 - quick fix

Until, the HMS cron server is able to run needed MySQL libs,
this is a hackish workaround.

Supervisor will be used to keep this script alive--at least for the weekend.
"""
from subprocess import Popen
import time

class DropboxPipelineWorkaround(object):

    def __init__(self, run_forever=True):

        if run_forever is True:
            while 1:
                self.run_workaround()
        else:
            self.run_workaround()

    def pause(self, num_minutes=10):
        print 'Pausing for %s minutes' % num_minutes
        num_seconds = num_minutes * 60
        time.sleep(num_seconds)

    def run_workaround(self):

        self.run_dropbox_retrieval_command()
        self.pause(10)

        self.run_pipeline_command()
        self.pause(10)

    def run_dropbox_retrieval_command(self):
        print 'Run Dropbox retrieval command'

        cmd_dropbox_retrieval = '/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/cron_scripts/get_dropbox_files_prod_hms.sh'
        cmd_args = cmd_dropbox_retrieval.split()

        p = Popen(cmd_args,
                shell=False,
                stdin=None,
                stdout=None,
                stderr=None,
                close_fds=True)


    def run_pipeline_command(self):
        print 'Run pipeline command'

        cmd_pipeline = '/www/gentb.hms.harvard.edu/code/gentb-site/gentb_website/cron_scripts/run_pipeline_prod_hms.sh'

        cmd_args = cmd_pipeline.split()

        p = Popen(cmd_args,
                shell=False,
                stdin=None,
                stdout=None,
                stderr=None,
                close_fds=True)



if __name__ == '__main__':
    temp_cmd_runner = DropboxPipelineWorkaround(run_forever=False)
