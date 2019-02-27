"""
Initialize a pipeline run from a pipeline and run it
"""

import time

from django.core.management.base import BaseCommand
from apps.pipeline.models import Pipeline


class PipelineDoesNotExist(Exception):
    pass


class Command(BaseCommand):
    help = __doc__

    def __print_status(self, program_run):
        pgruns = program_run.all_programs()
        for pgr in pgruns:
            print(pgr, pgr.run_time(), '\n')

    def add_arguments(self, parser):
            parser.add_argument(
                'pipeline_name',
                type=str,
                help='Name of the pipeline to test'
            )

            parser.add_argument(
                'update_rate',
                type=float,
                help='Rate to update status of program runs'
            )

    def handle(self, **kwargs):
        if 'update_rate' not in kwargs:
            update_rate = 1
        else:
            update_rate = kwargs['update_rate']

        pname = kwargs['pipeline_name']
        ps = Pipeline.objects.filter(name=pname)

        # TODO make exception more descriptive
        if not ps:
            raise PipelineDoesNotExist

        for pipeline in ps:
            prun = pipeline.run(for_test=True)

            # TODO see if this can update based on when program runs finish,
            # not on a clock
            # TODO more descriptive output
            while True:
                if prun.update_all():
                    self.__print_status(prun)

                    break

                self.__print_status(prun)

                # stop as soon as an error is found
                if prun.programs.filter(is_error=True).count() != 0:
                    break

                time.sleep(update_rate)
