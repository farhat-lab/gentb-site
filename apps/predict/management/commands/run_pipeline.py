
import sys
import time

from django.core.management.base import BaseCommand, CommandError

from apps.predict.models import PredictStrain
from apps.pipeline.models import ProgramRun

class Command(BaseCommand):
    help = """Schedule each of the pipeline tasks with the shell or LSF"""

    def handle(self, **options):
        for strain in PredictStrain.objects.filter(piperun__isnull=True):
            sys.stderr.write("Strain: %s, " % str(strain))
            dl = strain.check_download()
            if not dl:
                sys.stderr.write("NOT-READY [SKIP]\n")
                continue
            if strain.run():
                sys.stderr.write("SUBMITTED [OK]\n")
            else:
                sys.stderr.write("SUBMITTED [ERROR]\n")
            time.sleep(1)

        sys.stderr.write("\n")

        for strain in PredictStrain.objects.filter(
                piperun__programs__is_submitted=False,
                piperun__programs__is_complete=True,
                piperun__programs__is_error=True,
                ).distinct():
            sys.stderr.write("Re-Strain: %s " % str(strain))
            if strain.piperun.rerun():
                sys.stderr.write("SUBMITTED [OK]\n")
            else:
                sys.stderr.write("SUBMITTED [ERROR]\n")
            time.sleep(1)

        for piperun in ProgramRun.objects.filter(is_submitted=True, is_complete=False, is_error=False):
            sys.stderr.write("Run Status: %s " % str(piperun))
            if piperun.update_status():
                sys.stderr.write(" [COMPLETE]\n")
            elif piperun.is_error:
                sys.stderr.write(" [ERROR]\n")
            else:
                sys.stderr.write(" [WAITING]\n")

