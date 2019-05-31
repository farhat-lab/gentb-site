"""
Manage pipeline jobs, re-submit or submit as needed
"""
import sys

from datetime import datetime
from django.core.management.base import BaseCommand

from chore import JobSubmissionError

from apps.pipeline.models import PipelineRun, ProgramRun

def log(msg, *args, **kwargs):
    """Write consistantly to output"""
    kwargs['dt'] = datetime.now().isoformat()
    msg = "{dt}: " + msg
    sys.stderr.write(msg.format(*args, **kwargs) + "\n")

class Command(BaseCommand):
    """Schedule each of the pipeline tasks with the shell or LSF"""
    help = __doc__

    def handle(self, **options):
        """Called from the command line"""

        # These items were submitted but not complete yet, check status.
        for piperun in PipelineRun.objects.filter(
                programs__is_submitted=True,
                programs__is_complete=False,
                programs__is_error=False):

            piperun.update_all()

            log("STAT: {} |{}|".format(piperun, piperun.text_status()))

        # Rerun process
        for run in ProgramRun.objects.filter(is_submitted=False, is_error=True):
            if not run.has_input:
                continue
            log("RERUN: {}", run)
            run.error_text = ''
            run.is_error = False
            run.is_started = False
            run.is_complete = False
            run.is_submitted = True
            run.job_clean()
            try:
                run.job_submit(run.debug_text)
            except JobSubmissionError as err:
                run.is_error = True
                run.error_text = "Error while RE-SUBMITTING the job: " + str(err)
            except Exception as err: # pylint: disable=broad-except
                run.is_error = True
                run.error_text = "Owch Error: " + str(err)
            if run.error_text:
                log(run.error_text)
            run.save()
