#
# Copyright (C) 2017 Maha Farhat
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

import os
import time
import signal
import tempfile

from django.test import TestCase, override_settings
from autotest.base import ExtraTestCase

from apps.pipeline.models import Program, ProgramFile, Pipeline
from apps.pipeline.method import get_job_manager

DIR = os.path.dirname(__file__)
FIX = os.path.join(DIR, 'fixtures')

class JobManagerTest(TestCase):
    def setUp(self):
        self.manager = get_job_manager('apps.pipeline.method.shell')
        self.fn = tempfile.mktemp(prefix='test-job-')

    def tearDown(self):
        for x in range(20):
            fn = "%s.%d" % (self.fn, x)
            if os.path.isfile(fn):
                os.unlink(fn)
        if os.path.isfile(self.fn):
            os.unlink(self.fn)

    def test_shell_run(self):
        """Test that jobs can be run via the shell"""
        self.manager.submit('sleep_test_1', 'sleep 60')
        data = self.manager.status('sleep_test_1')
        self.assertIn(data['status'], ('sleeping', 'running'))
        self.manager.stop('sleep_test_1')
        data = self.manager.status('sleep_test_1')
        self.assertEqual(data['status'], 'stopped')

    def test_non_existant_id(self):
        """what happens when the job doesn't exist"""
        data = self.manager.status('sleep_test_0')
        self.assertEqual(data, {})
        self.manager.stop('sleep_test_0')

    def test_error_output(self):
        """When a command returns a non-zero status"""
        self.manager.submit('no_cmd', 'fidly dee')
        data = self.manager.status('no_cmd')
        while data['status'] == 'running':
            data = self.manager.status('no_cmd')
        self.assertEqual(data['status'], 'finished')
        self.assertEqual(data['error'], '/bin/sh: 1: fidly: not found')
        self.assertEqual(data['return'], 127)

    def assertDependantJobs(self, *cmds, **kw):
        for x, cmd in enumerate(cmds):
            cmd = 'sleep 0.1 && ' + (cmd % kw)
            depends = 'a%d' % (x - 1) if x else None
            self.manager.submit('a%d' % x, cmd, depends=depends)

        if 'call' in kw:
            kw['call']()

        ret = []
        job = 0
        timeout = len(cmds) * 10
        while job < len(cmds):
            data = self.manager.status('a%d' % job)
            if data.get('status', None) in ('finished', 'stopped', None):
                ret.append("%s:%s" % (data.get('status', 'No'), str(data.get('return', -1))))
                job += 1
                continue
            time.sleep(0.1)
            timeout -= 1
            self.assertTrue(timeout > 0, "Timeout waiting for dependant job %s" % cmds[job])

        return ret

    def test_dependant_jobs(self):
        """When one job needs a first job to complete"""
        ss = self.assertDependantJobs(
            'ls --help > %(fn)s.1',
            'grep OK %(fn)s.1 > %(fn)s.2',
            'wc %(fn)s.2 > %(fn)s.3',
            fn=self.fn)

        self.assertEqual(tuple(ss), ('finished:0',)*3)
        with open(self.fn+'.1', 'r') as fhl:
            self.assertTrue('Usage' in fhl.read())
        with open(self.fn+'.2', 'r') as fhl:
            self.assertEqual(fhl.read(), ' 0  if OK,\n')
        with open(self.fn+'.3', 'r') as fhl:
            self.assertEqual(fhl.read(), ' 1  3 11 %s.2\n' % self.fn)

    def test_dependant_error(self):
        """When the first job causes an error"""
        ss = self.assertDependantJobs(
            'ls %(fn)s.0 > %(fn)s.1',
            'ls %(fn)s.1 > %(fn)s.2',
            'ls %(fn)s.2 > %(fn)s.3',
            fn=self.fn)
        self.assertEqual(tuple(ss), ('finished:2', 'stopped:1', 'stopped:1'))

    def test_dependant_stopped(self):
        """When the first job is stoppped whole chain is stopped"""
        def stop():
            time.sleep(0.1)
            self.manager.stop('a0')
        ss = self.assertDependantJobs('sleep 60', 'sleep 1', 'sleep 1', call=stop)
        self.assertEqual(tuple(ss), ('stopped:9', 'stopped:1', 'stopped:1'))


class ProgramTest(ExtraTestCase):
    def setUp(self):
        super(ProgramTest, self).setUp()
        for name in ('one', 'two'):
            pf = ProgramFile.objects.create(name='file', store='test_%s.txt' % name)
            setattr(self, name, pf)
            setattr(self, name+'_fn', unicode(pf.store.file))

        self.program = Program(name='test', keep=False,
                command_line='ls -l ${file} > @{file}')
        self.program.save()
        self.program.files = [self.one]

    def test_no_input_error(self):
        """Test input error"""
        self.program.files = []
        with self.assertRaises(ValueError):
            dict(self.program.prepare_files())

    def test_no_output_error(self):
        """Test output error"""
        dict(self.program.prepare_files(file=self.one_fn))
        self.program.command_line = 'ls -l ${file} > @{output}'

        with self.assertRaises(ValueError):
            dict(self.program.prepare_files(file=self.one_fn))

    def test_io_output(self):
        """Test the processing of io"""
        self.program.command_line = 'A @{B} C ${D} E @{F} G ${H}'
        inputs = [('$', '', 'D', '', 9, 13), ('$', '', 'H', '', 23, 27)]
        outputs = [('@', '', 'B', '', 2, 6), ('@', '', 'F', '', 16, 20)]
        # Inputs are always first in the io
        self.assertEqual(list(self.program.io()), inputs + outputs)
        self.assertEqual(list(self.program.io(outputs=False)), inputs)
        self.assertEqual(list(self.program.io(inputs=False)), outputs)

    def test_program(self):
        """Test the program"""
        output = '/tmp/test_one.txt'
        files = dict(self.program.prepare_files(output_dir='/tmp'))
        self.assertEqual(files[('$', 'file', 6, 13)], self.one_fn)
        self.assertEqual(files[('@', 'file', 16, 23)], output)

        cmd = self.program.prepare_command(files)
        self.assertEqual(cmd, "ls -l %s > %s" % (self.one_fn, output))

    def test_prefix_suffix(self):
        """Test the use of prefix and suffix in command"""
        output = '/tmp/out_one.ls'
        self.one.name = 'foo'
        self.one.save()
        self.program.command_line = 'ls -l $test_{foo}.txt > @out_{foo}.ls'
        files = dict(self.program.prepare_files(output_dir='/tmp'))
        self.assertEqual(files[('@', 'foo', 24, 37)], output)
        cmd = self.program.prepare_command(files)
        self.assertEqual(cmd, "ls -l %s > %s" % (self.one_fn, output))

    def test_bin_directory(self):
        """Test the use of the custom BIN directory"""
        with self.settings(PIPELINE_BIN=FIX):
            self.program.command_line = 'sh ${bin}test.sh ${file}'
            for fn in dict(self.program.prepare_files(output_dir='/tmp')).values():
                self.assertTrue(os.path.isfile(fn), "File doesn't exist: %s" % fn)

    def test_brand_new_output(self):
        """Test the use of outputs with a new name"""
        self.program.command_line = 'ls > @{foo}.sam'
        files = dict(self.program.prepare_files(output_dir='/tmp', foo='gah.txt'))
        cmd = self.program.prepare_command(files)

    def test_command_combination(self):
        """New lines are considered command combinators"""
        self.program.command_line = 'ls\nwc\nls'
        files = dict(self.program.prepare_files(output_dir='/tmp'))
        cmd = self.program.prepare_command(files)
        self.assertEqual(cmd, "ls && wc && ls")

    def test_prepare_from_list(self):
        self.program.files = [self.one, self.two]
        self.program.command_line = 'ls ${file}_one.txt ${file}_two.txt '\
                                     + '${file}_one.ps ${file}_two.ps'
        files = list(self.program.prepare_files(output_dir='/tmp',
                                file=['/tmp/in_one.ps', '/tmp/in_two.ps']))

        self.assertEqual(files, [
            (('$', 'file', 3, 18), self.one_fn),
            (('$', 'file', 19, 34), self.two_fn),
            (('$', 'file', 35, 49), '/tmp/in_one.ps'),
            (('$', 'file', 50, 64), '/tmp/in_two.ps'),
        ])


class PipelineTest(ExtraTestCase):
    def setUp(self):
        self.pipeline = Pipeline.objects.create(name='TEST')
        self.file = ProgramFile(name='file', store='test_one.txt')
        self.fn = unicode(self.file.store.file)
        self.dir = os.path.dirname(self.fn)

    def setupPipeline(self, *programs):
        """Setup a list of programs in the pipeline"""
        wait = 'sleep 0.01 && '
        self.programs = [Program.objects.create(name=a, command_line=wait + b)
            for a, b in programs]
        for pos, program in enumerate(self.programs):
            self.pipeline.programs.create(order=pos, program=program)

    def assertProgram(self, result, inputs, outputs, data=None):
        def content(fn):
            with open(fn, 'r') as fhl:
                return fhl.read()
        outputs = outputs if isinstance(outputs, list) else [outputs]
        inputs = inputs if isinstance(inputs, list) else [inputs]
        self.assertEqual(result.input_files, "\n".join(inputs))
        self.assertEqual(result.output_files, "\n".join(outputs))
        if data is not None:
            ret = '---'.join([content(fn) for fn in outputs])
            self.assertEqual(ret, data)

        self.assertTrue(result.is_submitted)
        self.assertTrue(result.is_complete)
        self.assertEqual(result.error_text, None)
        self.assertFalse(result.is_error)

        self.assertGreater(result.input_size, 0)
        self.assertGreater(result.output_size, 0)
        self.assertGreater(result.duration, 0)

    @override_settings(PIPELINE_MODULE='apps.pipeline.method.fake')
    def test_pipeline(self):
        """Test the pipeline generation"""
        self.setupPipeline(
          ('A', 'ls -l ${file}.txt > @{file}.ls'),
          ('B', 'wc ${file}.ls > @{file}.c'),
          ('C', 'ls -l ${file}.txt ${file}.ls ${file}.c > @{file}.out'),
          ('D', 'wc ${file}.out > @{file}.out'),
        )
        result = self.pipeline.run("pipe", output_dir=self.dir, file=self.fn)
        limit = 40
        while not result.update_all():
            if limit == 40:
                get_job_manager().run_all()
            elif limit == 35:
                get_job_manager().finish_all()
            time.sleep(0.2)
            limit -= 1
            self.assertTrue(limit > 0, "Pipeline test timed out.")

        results = list(result.programs.all())
        self.assertEqual(len(results), 4)
        fn = self.fn[:-4] + ".%s"
        self.assertProgram(results[0], self.fn, fn % "ls")
        self.assertProgram(results[1], fn % "ls", fn % "c")
        self.assertProgram(results[2], [fn % "ls", self.fn, fn % "c"], fn % "out")
        self.assertProgram(results[3], fn % "out", fn % "out")

    @override_settings(PIPELINE_MODULE='apps.pipeline.method.fake')
    def test_duration(self):
        """Test the duration during a run"""
        self.setupPipeline(('DUR', 'sleep 5'))
        results = self.pipeline.run("pipe", output_dir=self.dir)
        result = results.programs.get()
        get_job_manager().run_all()
        for x in range(5):
            result.update_status()
            self.assertEqual(int(result.duration), x+1)
            time.sleep(1)


