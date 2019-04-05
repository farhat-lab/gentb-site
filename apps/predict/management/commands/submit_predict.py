"""
Submit the pipeline job as soon as the file download is complete.
"""
import sys
import time

from datetime import datetime
from django.core.management.base import BaseCommand

from chore import JobSubmissionError

from apps.pipeline.models import PrepareError
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

        try:
            return strain.run()
        except JobSubmissionError:
            raise ValueError("Can't run job")

    def handle(self, **_):
        """Called from the command line"""
        # Limit all interactions to a timeout (usually a few weeks)
        qset = PredictStrain.objects.filter(dataset__created__gt=get_timeout())

        for strain in qset.filter(piperun__isnull=True):
            try:
                self.submit_strain_pipeline(strain)
                log("RUN: {} ({})", strain, strain.pipeline)
                time.sleep(0.25)
            except PrepareError as err:
                log("BAD: {}: {}".format(strain, err))
            except IOError as err:
                log("ERR: {} (Bad Download) {}", strain, err)
            except Exception as err: # pylint: disable=broad-except
                log("ERR: {} ({}): {}", strain, strain.pipeline, err)
