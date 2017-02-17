
import sys
import time

from django.core.management.base import BaseCommand, CommandError

from apps.predict.models import PredictStrain

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

