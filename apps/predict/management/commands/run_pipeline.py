
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
            for input_file in (strain.file_one, strain.file_two):
                if input_file:
                    if input_file.retrieval_error:
                        sys.stderr.write("DOWNLOAD-ERROR [SKIP]\n")
                        continue
                    elif not input_file.retrieval_end:
                        sys.stderr.write("NOT-READY [SKIP]\n")
                        continue
            if strain.run():
                sys.stderr.write("SUBMITTED [OK]\n")
            else:
                sys.stderr.write("SUBMITTED [ERROR]\n")
            time.sleep(1)

        sys.stderr.write("\n")

        for piperun in ProgramRun.objects.filter(is_submitted=True, is_complete=False, is_error=False):
            sys.stderr.write("Run Status: %s " % str(piperun))
            if piperun.update_status():
                sys.stderr.write(" [COMPLETE]\n")
            elif piperun.is_error:
                sys.stderr.write(" [ERROR]\n")
            else:
                sys.stderr.write(" [WAITING]\n")

