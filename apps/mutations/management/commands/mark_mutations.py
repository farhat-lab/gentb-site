
from django.core.management.base import BaseCommand, CommandError

from apps.mutations.models import *

import sys
import logging
import argparse
LOGGER = logging.getLogger('apps.mutations')

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('infile', type=argparse.FileType('r'),
            default=sys.stdin, nargs='?',
            help='File with SNP names, one per line.')
        parser.add_argument('--unmark', action='store_true',
            help='Text files to load into the mutations app')

    def handle(self, infile, unmark=False, **kw):
        batch = []
        for snp in infile.readlines():
            snp = snp.strip() if snp else ''
            if not snp:
                break
            try:
                m = Mutation.objects.get(name=snp)
                m.predictor = not unmark
                m.save()
                sys.stderr.write('.')
                sys.stderr.flush()
            except Mutation.DoesNotExist:
                sys.stderr.write("%s Does not Exist\n" % snp)
        sys.stderr.write("\nDONE\n")


