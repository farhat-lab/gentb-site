# Generated by Django 2.2.13 on 2022-01-24 19:21

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('maptables', '0005_auto_20210913_1306'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mapdisplay',
            name='fill_max',
            field=models.DecimalField(decimal_places=4, default=-1.0, help_text='The maximum value this field will be, if set to -1 the value will be automaticly the maximum value in the set.', max_digits=20),
        ),
    ]
