
# pylint: disable=no-init,no-self-use
"""
Run a small pipeline on the configured shell.
"""

from django.core.management.base import BaseCommand, CommandError
from apps.pipeline.method import JobManager

class Command(BaseCommand):
    JOB_ID = 'sleep_60'

    def handle(self, **options):
        print "Job Manager: %s" % JobManager.name
        print "Starting job: sleep 60"

        JobManager.submit(self.JOB_ID, 'sleep 60')
        print "Job is running: %s" % str(JobManager.status(self.JOB_ID))

        print "Force Stop job..." 
        JobManager.stop(self.JOB_ID)

        print "Job is stopped: %s" % str(not JobManager.status(self.JOB_ID))


