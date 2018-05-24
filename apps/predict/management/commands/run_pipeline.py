"""
Cron job command for running each of the prediction pipelines as needed.
"""
import sys
import time

from django.core.management.base import BaseCommand

from apps.predict.models import PredictStrain, get_timeout
from apps.pipeline.models import ProgramRun

class Command(BaseCommand):
    help = """Schedule each of the pipeline tasks with the shell or LSF"""

    def handle(self, **options):
        """Called from the command line"""
        skipped, errors = 0, 0
        for strain in PredictStrain.objects.filter(piperun__isnull=True,\
                dataset__created__gt=get_timeout()):
            dl = strain.files_status
            sys.stderr.write("Strain: {}, Files: {}".format(strain, dl))
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
            if progrun.update_status():
                sys.stderr.write(" [COMPLETE]\n")
                progrun.piperun.predictstrain_set.first().update_status()
            elif progrun.is_error:
                sys.stderr.write(" [ERROR]\n")
            else:
                if 'fatal' in str(progrun).lower():
                    sys.stderr.write(" [FAILED]\n")
                    progrun.is_error = True
                    progrun.save()
                elif progrun.created < get_timeout():
                    sys.stderr.write(" [TIMEOUT]\n")
                    progrun.error_text = str(progrun.error_text) + '[JOB TIMEOUT]'
                    progrun.is_error = True
                    progrun.save()
                else:
                    sys.stderr.write(" [WAITING]\n")
