"""
Cron job command for running each of the prediction pipelines as needed.
"""
import sys
import time

from datetime import datetime
from django.core.management.base import BaseCommand

from apps.predict.models import (
    PredictStrain, get_timeout, STATUS_WAIT, STATUS_START, STATUS_ERROR
)

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

    @staticmethod
    def submit_strain_pipeline(strain):
        """Submit a strain pipleine to the job queue"""
        status = strain.files_status
        if status in (STATUS_WAIT, STATUS_START):
            # Waiting or busy doing file upload
            return False
        if status is STATUS_ERROR:
            raise IOError("Download Error")

        return strain.run()

    def handle(self, **options):
        """Called from the command line"""
        # Limit all interactions to a timeout (usually a few weeks)
        qset = PredictStrain.objects.filter(dataset__created__gt=get_timeout())

        for strain in qset.filter(
                piperun__programs__is_submitted=True,
                piperun__programs__is_complete=False,
                piperun__programs__is_error=False).distinct():
            # These items were submitted but not complete yet, check status.

            strain.update_status()

            stat = ''.join(['_> =!!?!'[bitset(
                run.is_submitted,
                run.is_complete,
                run.is_error)] for run in strain.piperun.programs.all()])

            log("STAT: {} |{}|".format(strain, stat))

        for strain in qset.filter(piperun__isnull=True):
            try:
                self.submit_strain_pipeline(strain)
                log("RUN: {} ({})", strain, strain.pipeline)
                time.sleep(0.25)
            except IOError:
                log("ERR: {} (Bad Download)", strain)
            except Exception as err:
                log("ERR: {} ({})", strain, strain.pipeline)
