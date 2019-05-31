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
"""
Pipeline is constructed from a pipeline pattern:

  Pipeline
    - Program
    - Program
    - Program

Each pipeline can then be run with PipelineFiles and outputside files
to produce a mirror record of a piepline being run:

  PipelineRun
    - ProgramRun
    - ProgramRun
    - ProgramRun

Testing is done by including PipelineFiles as testing files.
"""

import re
import os
from datetime import timedelta
import random

import logging
LOGGER = logging.getLogger('apps.pipeline')

from django.db.models import (
    Model, Q, PositiveIntegerField, FileField, SlugField, DateTimeField,
    BooleanField, CharField, ForeignKey, TextField, ManyToManyField,
)
from model_utils.models import TimeStampedModel

from chore import get_job_manager, tripplet, JobSubmissionError

from django.core.urlresolvers import reverse
from django.conf import settings
from django.utils.timezone import now
from django.utils.text import slugify
from django.utils.timezone import now

from .utils import file_as_inputs

class PrepareError(ValueError):
    """Error when preparing a command fails"""

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

    def get_absolute_url(self):
        return reverse('pipeline:detail', kwargs={'pk': self.pk})

    def run(self, name=None, commit=True, rerun=False, for_test=False, **kwargs):
        """
        Run this pipeline for this named identifier,
        the id should be unique.
        """

        # if name isn't given, use the existing pipeline name and append a
        # random number
        if name is None:
            name = self.name + '_' + str(random.randint(0, 9999))

        if commit:
            (runner, created) = self.runs.get_or_create(name=slugify(name))
            if not created and rerun:
                runner.programs.filter(is_error=True).update(is_submitted=False)
        else:
            runner = PipelineRun(name=slugify(name), pipeline=self)

        if for_test:
            kwargs.update(file_as_inputs(self.test_files))
        runner.run(commit=commit, for_test=for_test, **kwargs)
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
                        r'{(?P<literal>"?)(?P<name>[-\w]+)"?\}(?P<suffix>[^\s;|>]*)')

    name = CharField(max_length=128)
    description = TextField(null=True, blank=True,
            help_text='Describe the program and what it does in detail.')

    command_line = TextField(
            help_text='Write the command line using replacement syntax '
            'for inputs and outputs')
    keep = BooleanField(default=True,
            help_text="Should the output files be kept or deleted.")
    files = ManyToManyField('ProgramFile', related_name='use_in', blank=True,
            help_text="Files used when running this program")
    test_files = ManyToManyField('ProgramFile', related_name='test_in',
            blank=True, help_text="Files to test just this program.")
    memory = CharField(max_length=16, default='1000M', choices=[
        ('500M', '500MB'),
        ('1M', '1GB'),
        ('2G', '2GB'),
        ('5G', '5GB'),
        ('10G', '10GB'),
        ('20G', '20GB'),
        ('50G', '50GB'),
        ('100G', '100GB'),
    ], help_text="Amount of memory to request")
    limit = CharField(max_length=12, null=True, blank=True,
                      help_text="Specify amount of time to limit this process, "
                                "[days]-[hours]:[minutes]:[seconds]")
    threads = PositiveIntegerField(default=1)

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
        errors, files_out = [], {}

        # First take the files stored against this specific program.
        files_in = file_as_inputs(self.files)

        bins = getattr(settings, 'PIPELINE_BIN', None)
        if bins is not None:
            files_in['bin'] = [os.path.join(bins, d) for d in os.listdir(bins)]

        if output_dir is None:
            output_dir = '/tmp'

        # When running tests, we include inputs from the test files too.
        if for_test is True:
            file_as_inputs(self.test_files, files_in)

        # Override both files and test_files with inputs from the caller.
        # This is so programs can have reasonable defaults when running as
        # well as having either tests for the program or the whole pipeline.
        for name, value in inputs.items():
            if not isinstance(value, (tuple, list)):
                value = [value]
            for file_in in value:
                files_in[name].append(file_in)

        for (_io, prefix, literal, name, suffix, start, end) in self.io():
            try:
                if literal:
                    fname = name + suffix
                else:
                    args = (_io, prefix, name, suffix)
                    fname = self.prepare_file(files_in, files_out, *args)

                if '/' not in fname:
                    fname = os.path.join(output_dir, fname)
                yield ((_io, name, start, end), fname)
            except ValueError as err:
                errors.append(str(err))

        if errors:
            raise PrepareError("Error preparing command: \n * " + \
                    "\n * ".join(errors))

    def prepare_file(self, files_in, files_out, io, prefix, name, suffix):
        # Generic finder pattern
        pattern = "^%s(?P<name>[^/]*)%s$" % (prefix, suffix)
        # Get the file input from a list of available files.
        # It's reversed so the last added file is first.
        filenames = list(reversed(files_in[name]))
        for fname in filenames:
            ret = re.match(pattern, os.path.basename(fname))
            if ret:
                files_out[name] = ret.groupdict()['name']
                if io == '$':
                    return fname

            # For the pair-ended files, the file input is actually not a real
            # file input, but it needs to be made available to the output
            elif io == '@' and name not in files_out:
                files_out[name] = fname

        if io != '@':
            if not filenames:
                raise ValueError(self.ER3 % name)
            raise ValueError(self.ER1 % (name, pattern))

        # Populate / generate filename for output
        if name in files_out:
            # Make a new filename from the input's middle section
            # and the prefix and suffix named above.
            return ''.join([prefix, files_out[name], suffix])

        raise ValueError("Output '%s' unmatched from inputs." % name)


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
                raise PrepareError("Can't find file %s in %s" % (str(key), str(files)))
            filename = files[key]
            if ' ' in filename:
                filename = '"%s"' % filename
            # Replace the command line section that matches with out filename
            (start, end) = match.span()
            cmd = cmd[:start] + filename + cmd[end:]

        # Do the line replacer after the matching to preserve file placements.
        cmd = cmd.replace('\r', '').replace('\n\n', '\n')
        cmd = cmd.replace(';\n', '<NL>')
        cmd = cmd.replace('\n', ' && ')
        cmd = cmd.replace('<NL>', '\n')

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
        if pk is None:
            pk = 0
        return {
            'program': self.program,
            'job_id': "run_%d_%s" % (pk, slugify(self.program.name)),
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
        return "%s (%s)" % (self.name, os.path.basename(self.store.name))


class PipelineRun(TimeStampedModel):
    """
    A pipeline run is a representation of a single time when a given
    pipeline was submitted to the cluster. Each program attached to
    the pipeline creates a ProgramRun object attached to this PipelineRun
    each ProgramRun is a program submitted to the cluster with any returned
    values or printed strings attached.
    """
    name = SlugField(max_length=128, db_index=True)
    pipeline = ForeignKey(Pipeline, related_name='runs')
    run_as_test = PositiveIntegerField(null=True, blank=True,\
        help_text="Every (x) days, run this pipeline-run as a test.")
    clean_files = TextField(null=True, blank=True,
        help_text="List of extra files to remove when job is finished")

    def __str__(self):
        return "Pipeline Run %s" % self.name

    def text_status(self):
        """Generate a textual status of the run"""
        def bitset(*args):
            """Turn a set of booleans into a bit array and then into an integer"""
            return int(''.join(reversed([str(int(bool(arg))) for arg in args])), 2)

        return ''.join(['_> =!!?!'[bitset(
            progrun.is_submitted,
            progrun.is_complete,
            progrun.is_error
        )] for progrun in self.programs.all()])

    def get_absolute_url(self):
        """Return a link to this pipeline run for review"""
        return reverse('pipeline:run', kwargs={'pk': self.pk})

    def clean_filenames(self):
        """Return a list of fully qualified cleanable filenames"""
        if (self.clean_files or "").strip():
            return self.clean_files.split("\n")
        return []

    def clean_the_files(self):
        """Deletes any of the files marked for cleaning"""
        for fname in self.clean_filenames():
            try:
                os.unlink(fname)
            except (OSError, IOError):
                pass

    def run(self, commit=True, **kwargs):
        """Run this pipeline run (creates ProgramRun objects)"""
        runs = []
        if not commit:
            self.test_programs = []

        if 'clean_files' in kwargs:
            self.clean_files = '\n'.join(kwargs['clean_files'])
            if commit:
                self.save()

        for pipe in self.pipeline.programs.all():
            if commit:
                run, _ = ProgramRun.objects.get_or_create(piperun=self, **pipe.prepare(self.pk))
            else:
                run = ProgramRun(piperun=self, **pipe.prepare(self.pk))
                self.test_programs.append(run)
            runs.append(run)

        for prev, run, foll in tripplet(runs):
            if not run.is_submitted:
                if not run.submit(commit=commit, previous=prev, follower=foll, **kwargs):
                    return False
            else:
                data = get_job_manager().status(run.job_id, clean=False)
                if data.get('finished', None) and data.get('return', 1) != 1:
                    raise JobSubmissionError("Existing job already failed.")

            # Sort out the filenames for the next call in the chain
            for package, filename in run.program.prepare_files(**kwargs):
                name = package[1]
                if name in kwargs:
                    if isinstance(kwargs[name], list):
                        kwargs[name].append(filename)
                    else:
                        kwargs[name] = [kwargs[name], filename]
                else:
                    kwargs[name] = [filename]

        return True

    def all_programs(self):
        """Returns all the program runs with unrun pipelines appended"""
        ret = []
        runs = dict((p.program_id, p) for p in self.programs.all())
        for pipe in self.pipeline.programs.all():
            ret.append(runs.get(pipe.program_id, pipe))
        return ret

    def stop_all(self, msg='All Stopped'):
        """Forcefully stop all processes in this pipeline run"""
        return all([program.stop(msg=msg) for program in self.programs.all()])

    def update_all(self):
        """
        Update all pipeline project runs with their running status returns
        True if all programs are complete. False if any are still running.
        """
        qset = self.programs.filter(Q(is_submitted=False) | Q(is_complete=False))
        if all([program.update_status() for program in qset]):
            if qset.count():
                # Clean up step for all programs
                for program in self.programs.filter(program__keep=False):
                    program.delete_output_files()
                self.clean_the_files()
            return True
        return False

    def get_errors(self):
        """Return true if programs in the pipeline have errors"""
        self.update_all()
        qset = self.programs.filter(is_error=True)
        if qset.count() == 0:
            return None
        return '\n---\n'.join(qset.values_list('error_text', flat=True))

P_LOG = None

class ProgramRun(TimeStampedModel):
    piperun = ForeignKey(PipelineRun, related_name='programs')
    program = ForeignKey(Program, related_name='runs')
    job_id = SlugField(max_length=255, help_text="Name or ID of the job in the cloud runner")

    previous_id = SlugField(max_length=255, null=True, blank=True,\
        help_text="Name or ID of the previous job we depend on")
    follower_id = SlugField(max_length=255, null=True, blank=True,\
        help_text="Name or ID of the next job that depends on this")

    is_submitted = BooleanField(default=False)
    is_started = BooleanField(default=False)
    is_complete = BooleanField(default=False)
    is_error = BooleanField(default=False)

    input_size = PositiveIntegerField(null=True, blank=True,\
        help_text='Size in kilobytes of all input files.')
    output_size = PositiveIntegerField(null=True, blank=True,\
        help_text='Size in kilobytes of all output files.')
    duration = PositiveIntegerField(null=True, blank=True,\
        help_text='Number of seconds to run.')

    submitted = DateTimeField(null=True, blank=True)
    started = DateTimeField(null=True, blank=True)
    completed = DateTimeField(null=True, blank=True)

    input_files = TextField(null=True, blank=True)
    output_files = TextField(null=True, blank=True)
    debug_text = TextField("Command and Debug", null=True, blank=True)
    error_text = TextField("Error", null=True, blank=True)

    class Meta:
        ordering = ['created']

    def __str__(self):
        return self.job_id

    def get_absolute_url(self):
        return reverse('pipeline:job', kwargs={'pk': self.pk})

    def run_time(self):
        if self.duration:
            return str(timedelta(seconds=self.duration))
        if self.started:
            return now() - self.started
        return "-"

    def wait_time(self):
        if self.started:
            return str(self.started - self.submitted)
        if self.submitted:
            return now() - self.submitted
        return "-"

    def stop(self, msg='Stopped'):
        """Stop this program from running"""
        if self.is_submitted and not self.is_complete:
            ret = get_job_manager().stop(self.job_id)
            self.is_error = True
            self.is_complete = True
            self.error_text = msg
            self.save()
            return ret
        return True

    def submit(self, commit=True, previous=None, follower=None, **kwargs):
        """Submit the job and capture any errors"""
        try:
            cmd = self.prepare_command(commit=commit, **kwargs)

            job_kwargs = {}
            if previous is not None:
                job_kwargs['depend'] = self.previous_id = previous.job_id
            if follower is not None:
                job_kwargs['provide'] = self.follower_id = follower.job_id

            self.is_submitted = True
            self.submitted = now()

            if commit:
                job_kwargs['memory'] = str(self.program.memory)
                job_kwargs['threads'] = str(self.program.threads)
                if self.program.limit:
                    job_kwargs['limit'] = str(self.program.limit)
                job_kwargs['wckey'] = slugify(self.piperun.pipeline.name)
                self.job_submit(cmd, **job_kwargs)
                self.save()

        except (JobSubmissionError, PrepareError) as err:
            self.is_error = True
            self.error_text = str(err)
            self.kwargs = {}
            for key in kwargs:
                if key in ('output_dir', 'previous', 'follower'):
                    continue
                if isinstance(kwargs[key], (str, unicode)):
                    self.kwargs[key] = [kwargs[key]]
                elif isinstance(kwargs[key], (list, tuple)):
                    self.kwargs[key] = kwargs[key]
            if commit:
                self.save()
            return False
        return True

    def prepare_command(self, **kwargs):
        """Submit this job to the configured Job Manager"""
        # Save all the input and output files into database
        fsi = dict(self.program.prepare_files(**kwargs))
        self.input_files = '\n'.join([fsi[k] for k in fsi if k[0] == '$'])
        self.output_files = '\n'.join([fsi[k] for k in fsi if k[0] == '@'])

        # Remove any output files previously used.
        self.delete_output_files()

        cmd = self.program.prepare_command(dict(fsi))
        self.debug_text = cmd
        return cmd

    def job_submit(self, cmd, **kwargs):
        """Actually submit job to job_manager"""
        job_manager = get_job_manager()
        job_manager.submit(self.job_id, cmd, **kwargs)

    def job_clean(self):
        """Remove old command files"""
        job_manager = get_job_manager()
        job_manager.job_clean_fn(self.job_id, 'out')
        job_manager.job_clean_fn(self.job_id, 'err')

    def update_status(self, commit=True):
        """Take data from the job manager and populate the database"""
        job_manager = get_job_manager()
        if self.is_submitted and not self.is_complete:
            dur = None
            data = job_manager.status(self.job_id, clean=False)
            age = now() - self.submitted

            if not data and age > timedelta(hours=1):
                # This usually means the job is so old that it's gone from
                # the job manager queue and we have no further information about it
                self.is_complete = True
                self.is_error = True
                self.error_text = "Job Disapeared from Job Queue"
                return

            if data.get('status', 'notfound') in ('finished',):
                if data['finished'] and data['started']:
                    dur = data['finished'] - data['started']
                    self.duration = dur.total_seconds() + int(dur.microseconds > 0)
                self.completed = data['finished']
                self.is_complete = True
                self.is_error = data['return'] != 0
                if data['error']:
                    self.error_text = data['error'][:10240] # Limit errors to 10k
                self.input_size = self.update_size(*self.input_fn) / 1024.0
                self.output_size = self.update_size(*self.output_fn) / 1024.0

            if data.get('started', None) is not None:
                if not self.is_started:
                    self.is_started = True
                    self.started = data['started']
                # Save the duration so far
                dur = now() - data['started']
                # Round up any microseconds, useful for testing non-zero time
                self.duration = dur.total_seconds() + int(dur.microseconds > 0)

            if data and self.previous_id:
                for prev in ProgramRun.objects.filter(job_id=self.previous_id):
                    if prev.is_error:
                        job_manager.stop(self.job_id)
                        self.is_error = True

            if commit:
                self.save()

        # We're going to force an error out of hiding.
        if self.is_complete and self.error_text == 'None':
            (_, error) = job_manager.job_read(self.job_id, 'err')
            if error is not None:
                self.error_text = "Broken JobID error: " + error
            else:
                self.error_text = "Lost error for {}".format(self.job_id)
            self.is_error = True
            self.save()

        return self.is_complete

    @property
    def has_output(self):
        """Returns true if all the expected output files exist"""
        return len(self.output_fn) == len(self.output_filenames())

    @property
    def output_fn(self):
        """Returns a list of output filenames"""
        if self.output_files is not None:
            return [fn for fn in self.output_filenames() if os.path.isfile(fn)]
        return []

    @property
    def has_input(self):
        """Returns true if all the expected input files exist"""
        return len(self.input_fn) == len(self.input_filenames())

    @property
    def input_fn(self):
        """Returns a list of output filenames"""
        if self.input_files is not None:
            return [fn for fn in self.input_filenames() if os.path.isfile(fn)]
        return []

    def input_filenames(self):
        """Return a list of fully qualified input filenames"""
        if self.input_files is not None and self.input_files.strip():
            return self.input_files.replace('\r', '').split("\n")
        return []

    def output_filenames(self):
        """Return a list of fully qualified output filenames"""
        if self.output_files is not None and self.output_files.strip():
            return self.output_files.replace('\r', '').split("\n")
        return []

    def update_size(self, *files):
        """Takes a list of files as a string and returns the size in Kb"""
        return 1024 + sum([os.path.getsize(fn) for fn in files])

    def delete_output_files(self):
        """Deletes any of the output files"""
        for fname in self.output_fn:
            try:
                os.unlink(fname)
            except (OSError, IOError):
                pass
