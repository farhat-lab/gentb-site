#
# Copyright (C) 2017 - Dr. Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the 
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Submit a job with files to the pipeline system.
"""

import sys
import logging
LOGGER = logging.getLogger('apps.predict')

from django.core.management.base import BaseCommand, CommandError

from apps.predict.models import *
from apps.pipeline.models import *
from apps.predict.forms import *

class Command(BaseCommand):
    def add_arguments(self, parser):
        self.tp_choices = zip(*PredictDataset.FILE_TYPES)[0]
        super(Command, self).add_arguments(parser)
        parser.add_argument('--name', type=str, help='Name of this submission')
        parser.add_argument('--file', '-1', default=[], action='append',
            dest='files', help='Add files to this prediction run.')
        parser.add_argument('--file_two', '-2', default=[], action='append',
            dest='files_two', help='For pair ended fastq files, specify the '
            'second file after the first one IN ORDER.')
        parser.add_argument('--type', type=str, dest='tp',
            help='Type of prediction to run on these files.',
            choices=tuple(self.tp_choices))

        pipes = Pipeline.objects.values_list("pk", "name")
        pp = ", ".join("[%d=%s]" % p for p in pipes)
        parser.add_argument('--pipeline', type=int, dest='pipeline',
            help='Manually specify the pipeline to use.\n%s\n' % pp)

    def handle(self, name, files, tp, pipeline=None, **kw):
        FormClass = UploadForm.get_form(tp)
        print "Selected prediction: %s" % FormClass.doc_title

        post_dict = {
            'name': name,
        }

        form = FormClass(data, {})



