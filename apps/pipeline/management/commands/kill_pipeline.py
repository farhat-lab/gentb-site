"""
Kill one or many running PipelineRuns
"""

import sys

from django.core.management.base import BaseCommand
from apps.pipeline.models import PipelineRun


class Command(BaseCommand):
    help = __doc__

    verbosity = 1
    prun_name = ''

    def add_arguments(self, parser):
        by_name = parser.add_argument_group('Kill by PipelineRun name')
        killall = parser.add_argument_group('Kill all running PipelineRuns')

        by_name.add_argument(
            '--kill',
            type=str,
            help='Kill the specified pipeline',
        )

        killall.add_argument(
            '--killall',
            help='Kill all running PipelineRuns',
            action='store_true'
        )

    def handle(self, **kwargs):
        self.verbosity = kwargs.get('verbosity', 1)
        self.prun_name = kwargs.get('kill', '')

        if kwargs.get('kill') and kwargs.get('killall'):
            raise ValueError('Cannot use both --kill and --killall at the same time')

        prs = []

        if kwargs.get('kill'):
            prs = PipelineRun.objects.filter(name=self.prun_name)
        elif kwargs.get('killall'):
            prs = PipelineRun.objects.filter()

        if not prs:
            raise ValueError('Must specify either --kill or --killall')

        for pr in prs:
            pr.stop_all()

        sys.exit()
