
import sys
import time

from django.core.management.base import BaseCommand, CommandError

from apps.predict.models import PredictStrain
from apps.pipeline.models import ProgramRun

class Command(BaseCommand):
    help = """Schedule each of the pipeline tasks with the shell or LSF"""

    def handle(self, **options):
        skipped, errors = 0, 0
        for strain in PredictStrain.objects.filter(piperun__isnull=True):
            sys.stderr.write("Strain: %s, " % str(strain))
            dl = strain.files_status
            if dl is 0:
                skipped += 1
            elif dl is 1:
                skipped += 1
            elif dl is 2:
                errors += 1
            elif strain.run():
                sys.stderr.write("SUBMITTED [OK]\n")
            else:
                sys.stderr.write("SUBMITTED [ERROR]\n")
            time.sleep(1)

        if skipped or errors:
            sys.stderr.write("DOWNLOADS WARN: %d SKIPPED, %d ERRORS\n" % (skipped, errors))

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

        for progrun in ProgramRun.objects.filter(is_submitted=True, is_complete=False, is_error=False):
            sys.stderr.write("Run Status: %s " % str(progrun))
            if 'fatal' in str(progrun).lower():
                progrun.is_error = True
                progrun.save()
            if progrun.update_status():
                sys.stderr.write(" [COMPLETE]\n")
                progrun.piperun.predictstrain_set.first().update_status()
            elif progrun.is_error:
                sys.stderr.write(" [ERROR]\n")
            else:
                sys.stderr.write(" [WAITING]\n")

