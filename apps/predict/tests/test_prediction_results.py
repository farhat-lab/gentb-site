#
# Copyright (C) 2021 Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Test parsing various prediction json file formats.
"""
import os
import json
from decimal import Decimal

from extratest.base import ExtraTestCase

from apps.pipeline.models import PipelineRun

from ..models import PredictDataset, PredictStrain
from ..predict_data import decypher_predict_format

class TestPredictFiles(ExtraTestCase):
    fixtures = ['drugs',]

    def setUp(self):
        self.dts = PredictDataset.objects.create(
            pk=101,
            md5='F00BA5',
            title='Testing Dataset',
            has_notified=True,
            file_directory=os.path.join(self.source_dir, 'predictions'),
        )

    def _piperun(self, filename):
        fn = os.path.basename(filename)
        piperun = PipelineRun.objects.create(
            name=f'Pipelin for Testing {fn}',
        )
        piperun.programs.create(
            job_id=f'progrun_{fn}',
            is_submitted=True,
            is_started=True,
            is_complete=True,
            output_files=f"{filename}",
            job_state='COMPLETED',
        )
        return piperun

    def _strain(self, filename):
        """Creates a strain from the given json file"""
        fn = os.path.basename(filename)
        return PredictStrain.objects.create(
            name=f'Strain for Testing {fn}',
            dataset=self.dts,
            pipeline=None,
            piperun=self._piperun(filename),
        )

    def get_file(self, *fn):
        filepath = os.path.join(self.source_dir, *fn)
        if not os.path.isfile(filepath):
            raise IOError(f"Can't find requested data file {filepath}")
        return filepath

    def test_piperun(self):
        fn = self.get_file('default', 'matrix.json')
        piperun = self._piperun(fn)
        program = piperun.programs.get()
        self.assertTrue(program.has_output)
        self.assertEqual(program.output_fn, [fn])

    def test_strain(self):
        fn = self.get_file('default', 'matrix.json')
        strain = self._strain(fn)
        self.assertEqual(list(strain.output_files), [(None, fn, 'matrix.json', -1)])
        self.assertTrue(strain.has_prediction)

    def test_raw_prediction(self):
        strain = self._strain(self.get_file('default', 'matrix.json'))
        for err, dat in strain.get_raw_prediction():
            self.assertNotEqual(err, None)
            metadata, data = decypher_predict_format(dat)
            self.assertEqual(metadata, {'dr': '0', 'drug_code': 'inh', 'fneg': 4.95, 'fpos': 8.76})
            self.assertEqual(data, [
                (None, None, None, None, None),
                (None, None, None, None, None),
                ('DEL_CF_2155296_d816G_273_katG', None, None, None, None),
                (None, None, None, None, None)
            ])
            break

    def test_generated_results(self):
        strain = self._strain(self.get_file('default', 'matrix.json'))
        self.assertEqual(strain.results.count(), 0)
        strain.generate_results()
        self.assertEqual(strain.results.count(), 13)
        self.assertEqual(strain.results.get(drug__code='INH').probability, 0.0)
        self.assertEqual(strain.results.get(drug__code='INH').false_positive, Decimal('8.76'))
        self.assertEqual(strain.results.get(drug__code='INH').false_negative, Decimal('4.95'))

