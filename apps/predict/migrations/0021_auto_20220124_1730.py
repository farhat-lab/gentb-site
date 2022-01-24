# Generated by Django 2.2.13 on 2022-01-24 22:30

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0020_auto_20220124_1555'),
    ]

    operations = [
        migrations.AddField(
            model_name='predictdataset',
            name='strains_count',
            field=models.IntegerField(default=0, verbose_name='Number of strains total'),
        ),
        migrations.AddField(
            model_name='predictdataset',
            name='strains_ready',
            field=models.IntegerField(default=0, verbose_name='Number of strains ready'),
        ),
        migrations.AlterField(
            model_name='predictdataset',
            name='status',
            field=models.CharField(choices=[('FILE_NONE', 'No Files Uploaded'), ('FILE_WAIT', 'Dataset Confirmed'), ('FILE_START', 'File Retrieval Started'), ('FILE_ERROR', 'File Retrieval Failed'), ('RUN_NONE', 'No Strains to process'), ('RUN_WAIT', 'File Retrieval Success'), ('RUN_START', 'Processing Started'), ('RUN_ERROR', 'Processing Failed'), ('RUN_DONE', 'Processing Success'), ('READY', 'Prediction Ready'), ('INVALID', 'Lacks Quality'), ('TIMEOUT', 'Processing Timed Out'), ('', 'Status Unknown')], default='', max_length=10, verbose_name='Status'),
        ),
    ]
