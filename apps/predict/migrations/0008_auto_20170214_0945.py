# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('uploads', '0001_initial'),
        ('pipeline', '0001_initial'),
        ('predict', '0007_auto_20170202_1539'),
    ]

    operations = [
        migrations.CreateModel(
            name='PredictPipeline',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('file_type', models.CharField(max_length=25, choices=[(b'vcf', b'Variant Call Format (VCF)'), (b'fastq', b'FastQ Single Ended Nucleotide Sequence'), (b'fastq-pair', b'FastQ Pair Ended Nucleotide Sequences'), (b'manual', b'Mutations Manual Entry')])),
                ('is_default', models.BooleanField(default=False)),
                ('pipeline', models.ForeignKey(to='pipeline.Pipeline')),
            ],
        ),
        migrations.CreateModel(
            name='PredictStrain',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=128)),
            ],
        ),
        migrations.RemoveField(
            model_name='datasetscriptrun',
            name='dataset',
        ),
        migrations.RemoveField(
            model_name='predictdataset',
            name='fastq_type',
        ),
        migrations.RemoveField(
            model_name='predictdataset',
            name='has_prediction',
        ),
        migrations.AlterField(
            model_name='predictdataset',
            name='file_type',
            field=models.CharField(max_length=25, choices=[(b'vcf', b'Variant Call Format (VCF)'), (b'fastq', b'FastQ Single Ended Nucleotide Sequence'), (b'fastq-pair', b'FastQ Pair Ended Nucleotide Sequences'), (b'manual', b'Mutations Manual Entry')]),
        ),
        migrations.DeleteModel(
            name='DatasetScriptRun',
        ),
        migrations.AddField(
            model_name='predictstrain',
            name='dataset',
            field=models.ForeignKey(related_name='strains', to='predict.PredictDataset'),
        ),
        migrations.AddField(
            model_name='predictstrain',
            name='file_one',
            field=models.ForeignKey(related_name='link_a', blank=True, to='uploads.UploadFile', null=True),
        ),
        migrations.AddField(
            model_name='predictstrain',
            name='file_two',
            field=models.ForeignKey(related_name='link_b', blank=True, to='uploads.UploadFile', null=True),
        ),
        migrations.AddField(
            model_name='predictstrain',
            name='pipeline',
            field=models.ForeignKey(to='pipeline.Pipeline'),
        ),
        migrations.AddField(
            model_name='predictstrain',
            name='piperun',
            field=models.ForeignKey(blank=True, to='pipeline.PipelineRun', null=True),
        ),
    ]
