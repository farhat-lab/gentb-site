"""
Take the csv/tsv data from coding and non-coding files and insert into
Gene Locus database.
"""

import logging
from collections import OrderedDict

from django.core.management.base import BaseCommand

LOGGER = logging.getLogger('apps.mutations')

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--files', nargs='+', type=int,
            help='Specify the code files')

    def handle(self, dataset_id=None, retry=False, **options):
        for f in files:
            if not os.path.isfile(f):
                raise IOError("File doesn't exist: %s" % f)

        for f in files:
            self.process_file(f)

    def process_file(self, filename):
        with open(filename, 'r') as fhl:
            self.process_data(fhl.read().split('\n'))

    @staticmethod
    def process_data(data):
        drugs = OrderedDict()
        drug = None
        for item in data:
            item = item.strip()
            if not item:
                continue
            if item and not ' ' in item:
                drug = item
                DRUGS[drug] = []
                continue
            if not drug:
                logging.warning("Drug missing in mutation list (%s).", item)
                continue

            (x, result) = item.split(' ', 1)

            drugs[drug].append(result)

