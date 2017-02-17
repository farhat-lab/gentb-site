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

import re

from os.path import getsize, isfile, basename, join
from collections import defaultdict
from datetime import datetime

from django.db.models import *
from model_utils.models import TimeStampedModel

from django.conf import settings
from django.utils.timezone import now
from django.utils.text import slugify

import logging
LOGGER = logging.getLogger('apps.pipeline')

from apps.pipeline.method import JobManager

class Pipeline(Model):
    """
    Keeps a list of programs to run and what they do.
    """
    name        = CharField(max_length=128)
    description = TextField(null=True, blank=True,
            help_text='Describe the pipeline and what it does in detail.')
    test_files  = ManyToManyField('ProgramFile', related_name='tested_in',
            help_text='Input files used to run the pipeline test', blank=True)

    def __str__(self):
        return self.name

    def run(self, name, **kwargs):
        """
        Run this pipeline for this named identifier,
        the id should be unique.
        """
        runner = self.runs.create(name=name)
        runner.run(**kwargs)
        return runner


class Program(Model):
    """
    A program is a single step in a pipeline, a program that's run.

    Example Command:

        command_line: samtools view -bS ${sam}.sam | samtools sort | samtools index ${sam}
        inputs: ${file}.sam
        outputs: @{file}.bam
          The {file} part in the original filename is taken from the input and used
          For the output filename. The output file may end up in a different directory.

        stampy.py -g ${fasta} -h ${fasta} -o @{file}.sam -f sam -M ${fastq1}.fastq ${fastq2}.fastq;

    """
    ER1 = "Input '%s' no file matches pattern '%s'"
    ER2 = "Input '%s' file '%s' doesn't exist."
    ER3 = "Input '%s' not available in inputs"

    PARSER = re.compile(r'(?P<io>\$|\@)(?P<prefix>[^{]*)' \
                        r'{(?P<name>\w+)\}(?P<suffix>[^\s;|>]*)')

    name = CharField(max_length=128)
    description = TextField(null=True, blank=True,
            help_text='Describe the program and what it does in detail.')
    requirements = TextField(null=True, blank=True,
            help_text='List of requirements, one per line')

    command_line = TextField(
            help_text='Write the command line using replacement syntax '
            'for inputs and outputs')
    keep = BooleanField(default=True,
            help_text="Should the output files be kept or deleted.")
    files = ManyToManyField('ProgramFile', related_name='use_in', blank=True,
            help_text="Files used when running this program")
    test_files = ManyToManyField('ProgramFile', related_name='test_in',
            blank=True, help_text="Files to test just this program.")

    def __str__(self):
        return self.name

    def io(self, outputs=True, inputs=True):
        """Parse the command line for input and output variables
        inputs ALWAYS come out first and then outputs.

          inputs  - Should inputs be yielded, set to false to ignore
          outputs - Should outputs be yielded, set to false to ignore

        Yields regex match items items.
        """
        # We will stash the outputs for the post-section so the
        outs = []
        for match in self.PARSER.finditer(self.command_line):
            data = match.groups() + match.span()
            if data[0] == '@' and outputs is True:
                outs.append(data)

            elif data[0] == '$' and inputs is True:
                yield data

        for data in outs:
            yield data

    def as_inputs(self, m2m, save_to=None):
        """Adds each m2m to the save_to dictionary"""
        if not save_to:
            save_to = defaultdict(list)
        for pf in m2m.all():
            save_to[pf.name].append(unicode(pf.store.file))
        return save_to

    def prepare_files(self, for_test=False, output_dir=None, **inputs):
        """
        Prepares the files by testing the command line's required filenames
        while causing an error if there are pipeline mismatches.
        
        This is the main call for making a ready to use command line from
        this program given a set of input filenames.

         - for_test   - If True, will use test_files as a source of input data
         - output_dir - Should be set to the location where output files need
                        to be saved to.
         - **inputs   - A dictionary of inputs named in the command line which
                        override any program files or test_files.

          Yields a (name, filename) tuple suitable for a dictionary of both
          input and output filenames.
        """
        errors, pout = [], {}

        # First take the files stored against this specific program.
        files_in = self.as_inputs(self.files)

        if output_dir is None:
            output_dir = '/tmp'

        # When running tests, we include inputs from the test files too.
        if for_test is True:
            self.as_inputs(self.test_files, files_in)

        # Override both files and test_files with inputs from the caller.
        # This is so programs can have reasonable defaults when running as
        # well as having either tests for the program or the whole pipeline.
        for name, value in inputs.items():
            if not isinstance(value, (tuple, list)):
                value = [value]
            for file_in in value:
                files_in[name].append(file_in)

        for (io, prefix, name, suffix, start, end) in self.io():
            if io == '$':
                # Get the file input from a list of available files.
                # It's reversed so the last added file is first.
                fns = list(reversed(files_in[name]))

                pattern = "%s(?P<name>[^/]+)%s" % (prefix, suffix)

                for fn in fns:
                    ret = re.match(pattern, basename(fn))
                    if ret:
                        yield ((io, name, start, end), fn)
                        pout[name] = ret.groupdict()['name']

                if not fns:
                    errors.append(self.ER3 % name)

                elif name not in pout:
                    print fns
                    errors.append(self.ER1 % (name, pattern))

            elif io == '@':
                # Populate / generate filename for output
                if name in pout:
                    # Make a new filename from the input's middle section
                    # and the prefix and suffix named above.
                    new_name = ''.join([prefix, pout[name], suffix])
                    fn = join(output_dir, new_name)
                    yield ((io, name, start, end), fn)
                else:
                    errors.append("Output '%s' unmatched from inputs." % name)

        if errors:
            raise ValueError("Error preparing command: \n * " + \
                    "\n * ".join(errors))

    def prepare_command(self, files):
        """
        Prepares the command line itself with the given files, all files must
        already be full paths prepared by prepare_files or a suitable test.
        """
        # Now we have inputs and outputs prepared, we can construct a command
        # With all the input and output filenames in the right places.
        cmd = self.command_line

        for match in reversed(list(self.PARSER.finditer(cmd))):
            data = match.groupdict()
            key = (data['io'], data['name']) + match.span()
            if key not in files:
                raise KeyError("Can't find file %s in %s" % (str(key), str(files)))
            filename = files[key]
            if ' ' in filename:
                filename = '"%s"' % filename
            # Replace the command line section that matches with out filename
            (start, end) = match.span()
            cmd = cmd[:start] + filename + cmd[end:]

        return cmd


class PipelineProgram(Model):
    pipeline = ForeignKey(Pipeline, related_name='programs')
    program = ForeignKey(Program, related_name='pipelines')
    order = PositiveIntegerField(default=0)

    class Meta:
        ordering = ['order']

    def __str__(self):
        return "Program '%s' for Pipeline '%s'" % (str(self.program), str(self.pipeline))

    def prepare(self, pk):
        """Get the pipeline program ready for running"""
        return {
            'program': self.program,
            'job_id': "run_%d_%s" % (pk, slugify(self.program.name))
        }


class ProgramFile(Model):
    """
    Any useful data file used by the pipeline.
    """
    name = SlugField(max_length=32)
    store = FileField('File', upload_to='pipeline/files')
    description = TextField(null=True, blank=True,
            help_text='Describe the file is and what it does in detail.')

    def __str__(self):
        return self.name


class PipelineRun(TimeStampedModel):
    name = SlugField(max_length=128, db_index=True)
    pipeline = ForeignKey(Pipeline, related_name='runs')

    def __str__(self):
        return "Pipeline Run %s" % self.name

    def run(self, **kwargs):
        for pipe in self.pipeline.programs.all():
            run = ProgramRun.objects.create(piperun=self, **pipe.prepare(self.pk))
            try:
                for (name, fn) in run.submit(**kwargs):
                    if name in kwargs:
                        if isinstance(kwargs[name], list):
                            kwargs[name].append(fn)
                        else:
                            kwargs[name] = [kwargs[name], fn]
                    else:
                        kwargs[name] = [fn]
            except ValueError as err:
                run.is_error = True
                run.error_text = str(err)
                run.save()
                return False
            kwargs['previous'] = run
        return True

    def update_all(self):
        """
        Update all pipeline project runs with their running status returns
        True if all programs are complete. False if any are still running.
        """
        return all([program.update_status()
            for program in self.programs.all()])


class ProgramRun(TimeStampedModel):
    piperun = ForeignKey(PipelineRun, related_name='programs')
    program = ForeignKey(Program, related_name='runs')
    job_id  = SlugField(max_length=255, help_text="Name or ID of the job in the cloud runner")
    previous_id = SlugField(max_length=255, help_text="Name or ID of the previous job we depend on")
    
    is_submitted = BooleanField(default=False)
    is_complete = BooleanField(default=False)
    is_error = BooleanField(default=False)

    input_size = PositiveIntegerField(null=True, blank=True,
            help_text='Size in kilobytes of all input files.')
    output_size = PositiveIntegerField(null=True, blank=True,
            help_text='Size in kilobytes of all output files.')
    duration = PositiveIntegerField(null=True, blank=True,
            help_text='Number of seconds to run.')

    input_files = TextField(null=True, blank=True)
    output_files = TextField(null=True, blank=True)
    error_text = TextField(null=True, blank=True)

    def __str__(self):
        return self.job_id

    def submit(self, previous=None, **kwargs):
        files = dict(self.program.prepare_files(**kwargs))
        
        # Save all the input and output files into database
        fs = files.items()
        self.input_files = '\n'.join([fn for (k, fn) in fs if k[0] == '$'])
        self.output_files = '\n'.join([fn for (k, fn) in fs if k[0] == '@'])

        if previous is not None:
            self.previous_id = previous.job_id

        cmd = self.program.prepare_command(files)
        if JobManager.submit(self.job_id, cmd, depends=self.previous_id):
            self.is_submitted = True
        else:
            raise ValueError("Job could not be submitted to Job Manager.")

        self.save()
        # Return processed files so outputs can be used
        # as inputs to the next command in the pipeline
        return [(a[1], b) for (a,b) in files.items()]

    def update_status(self):
        """Take data from the job manager and populate the database"""
        if self.is_submitted and not self.is_complete:
            data = JobManager.status(self.job_id, clean=True)
            if data and data['return'] is not None:
                # The duration is
                dur = data['finished'] - data['started']

                self.is_complete = True
                self.is_error = data['return'] != 0
                self.error_text = data['error']
                self.input_size = self.update_size(self.input_files) / 1024.0
                self.output_size = self.update_size(self.output_files) / 1024.0
            elif data:
                # Save the duration so far
                dur = datetime.now() - data['started']

            if data:
                # Round up any microseconds, useful for testing non-zero time
                self.duration = dur.total_seconds() + int(dur.microseconds > 0)
                self.save()

        return self.is_complete

    def update_size(self, files):
        """Takes a list of files as a string and returns the size in Kb"""
        return 1024 + sum([getsize(fn)
            for fn in files.split("\n") if isfile(fn)])

