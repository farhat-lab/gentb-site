# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import django.utils.timezone
import model_utils.fields


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Pipeline',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=128)),
                ('description', models.TextField(help_text=b'Describe the pipeline and what it does in detail.', null=True, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='PipelineProgram',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('order', models.PositiveIntegerField(default=0)),
                ('pipeline', models.ForeignKey(related_name='programs', to='pipeline.Pipeline')),
            ],
            options={
                'ordering': ['order'],
            },
        ),
        migrations.CreateModel(
            name='PipelineRun',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('name', models.SlugField(max_length=128)),
                ('pipeline', models.ForeignKey(related_name='runs', to='pipeline.Pipeline')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Program',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=128)),
                ('description', models.TextField(help_text=b'Describe the program and what it does in detail.', null=True, blank=True)),
                ('requirements', models.TextField(help_text=b'List of requirements, one per line', null=True, blank=True)),
                ('command_line', models.TextField(help_text=b'Write the command line using replacement syntax for inputs and outputs')),
                ('keep', models.BooleanField(default=True, help_text=b'Should the output files be kept or deleted.')),
            ],
        ),
        migrations.CreateModel(
            name='ProgramFile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.SlugField(max_length=32)),
                ('store', models.FileField(upload_to=b'pipeline/files', verbose_name=b'File')),
                ('description', models.TextField(help_text=b'Describe the file is and what it does in detail.', null=True, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='ProgramRun',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, verbose_name='created', editable=False)),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, verbose_name='modified', editable=False)),
                ('job_id', models.SlugField(help_text=b'Name or ID of the job in the cloud runner', max_length=255)),
                ('previous_id', models.SlugField(help_text=b'Name or ID of the previous job we depend on', max_length=255)),
                ('is_submitted', models.BooleanField(default=False)),
                ('is_complete', models.BooleanField(default=False)),
                ('is_error', models.BooleanField(default=False)),
                ('input_size', models.PositiveIntegerField(help_text=b'Size in kilobytes of all input files.', null=True, blank=True)),
                ('output_size', models.PositiveIntegerField(help_text=b'Size in kilobytes of all output files.', null=True, blank=True)),
                ('duration', models.PositiveIntegerField(help_text=b'Number of seconds to run.', null=True, blank=True)),
                ('input_files', models.TextField(null=True, blank=True)),
                ('output_files', models.TextField(null=True, blank=True)),
                ('error_text', models.TextField(null=True, blank=True)),
                ('piperun', models.ForeignKey(related_name='programs', to='pipeline.PipelineRun')),
                ('program', models.ForeignKey(related_name='runs', to='pipeline.Program')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='program',
            name='files',
            field=models.ManyToManyField(help_text=b'Files used when running this program', related_name='use_in', to='pipeline.ProgramFile', blank=True),
        ),
        migrations.AddField(
            model_name='program',
            name='test_files',
            field=models.ManyToManyField(help_text=b'Files to test just this program.', related_name='test_in', to='pipeline.ProgramFile', blank=True),
        ),
        migrations.AddField(
            model_name='pipelineprogram',
            name='program',
            field=models.ForeignKey(related_name='pipelines', to='pipeline.Program'),
        ),
        migrations.AddField(
            model_name='pipeline',
            name='test_files',
            field=models.ManyToManyField(help_text=b'Input files used to run the pipeline test', related_name='tested_in', to='pipeline.ProgramFile', blank=True),
        ),
    ]
