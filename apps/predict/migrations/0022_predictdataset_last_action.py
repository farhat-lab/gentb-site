# Generated by Django 2.2.13 on 2022-01-24 23:03

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('predict', '0021_auto_20220124_1730'),
    ]

    operations = [
        migrations.AddField(
            model_name='predictdataset',
            name='last_action',
            field=models.DateTimeField(blank=True, null=True, verbose_name='Last Completed Action'),
        ),
    ]
