"""
Manage pipeline jobs, re-submit or submit as needed
"""
import sys

from datetime import datetime
from django.core.management.base import BaseCommand

from chore import JobSubmissionError

from apps.pipeline.models import ProgramRun

def bitset(*args):
    """Turn a set of booleans into a bit array and then into an integer"""
    return int(''.join(reversed([str(int(bool(arg))) for arg in args])), 2)

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

        # Rerun process
        for run in ProgramRun.objects.filter(is_submitted=False, is_error=True):
            if not run.has_input:
                continue
            log("RERUN: {}", run)
            run.is_error = False
            run.is_started = False
            run.is_finished = False
            run.is_submitted = True
            try:
                run.job_submit(run.debug_text)
            except JobSubmissionError as err:
                run.is_error = True
                run.error_text = "Error while RE-SUBMITTING the job: " + str(err)
            except Exception as err:
                run.is_error = True
                run.error_text = "Owch Error: " + str(err)
            run.save()
