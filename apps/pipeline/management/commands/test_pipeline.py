"""
Initialize a pipeline run from a pipeline and run it
"""

import sys
import time
import datetime

from django.core.management.base import BaseCommand
from apps.pipeline.models import Pipeline


class PipelineDoesNotExist(Exception):
    pass


class Command(BaseCommand):
    help = __doc__

    verbosity = 1
    pipeline_name = ''
    pipeline_run_name = ''

    def __print_status(self, program_run):
        labels = ('PROGRAM', 'TIME', 'TOTAL')

        # data for each row
        rows = []

        # longest string in each column
        cmax = {}

        for l in labels:
            cmax[l] = len(l)

        prev_seconds = 0

        pgruns = program_run.all_programs()
        for pgr in pgruns:
            pgr_name = str(pgr)
            pgr_total = pgr.run_time()

            # convert time to seconds
            pgr_seconds = 0
            if pgr_total != '-':
                parts = map(int, pgr_total.split(':'))
                parts = map(lambda a, b: a * b, parts, [3600, 60, 1])
                pgr_seconds = reduce(lambda a, b: a + b, parts, pgr_seconds)

            # find delta from previous program's time
            pgr_delta = '-'
            if pgr_seconds > prev_seconds or (pgr.completed and not pgr.is_error):
                pgr_delta = str(datetime.timedelta(seconds=pgr_seconds - prev_seconds))

            prev_seconds = pgr_seconds

            # max column lengths
            cmax['PROGRAM'] = max(cmax['PROGRAM'], len(pgr_name))
            cmax['TIME'] = max(cmax['TIME'], len(pgr_delta))
            cmax['TOTAL'] = max(cmax['TOTAL'], len(pgr_total))

            rows.append((pgr_name, pgr_delta, pgr_total))

        # template widths defined by longest output per column
        template = ' '.join(['{' + str(i) + ':%d}' for i in range(len(labels))])
        template = template % tuple(cmax[l] + 1 for l in labels)

        # print it all out
        # clears screen and puts cursor at top left
        if self.verbosity < 2:
            sys.stdout.write('\033[2J\033[1;1H')

        print('--- {} ---\n'.format(self.pipeline_run_name))
        print(template.format(*labels))
        for row in rows:
            print(template.format(*row))
        print('')

    def add_arguments(self, parser):
        parser.add_argument(
            'pipeline_name',
            type=str,
            help='Name of the pipeline to test'
        )

        parser.add_argument(
            '--update_rate',
            type=float,
            help='Rate to update status of program runs',
            default=1
        )

    def handle(self, **kwargs):
        self.pipeline_name = kwargs.get('pipeline_name')
        self.verbosity = kwargs.get('verbosity', 1)

        update_rate = kwargs.get('update_rate')
        pname = self.pipeline_name

        ps = Pipeline.objects.filter(name=pname)

        # TODO make exception more descriptive
        if not ps:
            raise PipelineDoesNotExist

        # TODO this should probably only handle one pipeline
        for pipeline in ps:
            prun = pipeline.run(for_test=True)
            self.pipeline_run_name = prun.name

            # TODO see if this can update based on when program runs finish,
            # not on a clock
            while True:
                if prun.update_all():
                    self.__print_status(prun)

                    break

                self.__print_status(prun)

                # stop as soon as an error is found
                if prun.programs.filter(is_error=True).count() != 0:
                    break

                time.sleep(update_rate)
