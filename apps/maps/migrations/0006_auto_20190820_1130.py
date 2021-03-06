# Generated by Django 2.0.13 on 2019-08-20 15:30

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('maps', '0005_countryhealth'),
    ]

    operations = [
        migrations.AlterField(
            model_name='country',
            name='region',
            field=models.IntegerField(blank=True, choices=[(0, 'Other'), (1, 'World'), (2, 'Agrica'), (9, 'Oceania'), (19, 'Americas'), (21, 'North America'), (142, 'Asia'), (150, 'Europe'), (419, 'Latin America and the Caribbean')], null=True),
        ),
        migrations.AlterField(
            model_name='country',
            name='subregion',
            field=models.IntegerField(blank=True, choices=[(14, 'Eastern Africa'), (17, 'Middle Africa'), (15, 'Northern Africa'), (18, 'Southern Africa'), (11, 'Western Africa'), (29, 'Caribbean'), (13, 'Central America'), (5, 'South America'), (143, 'Central Asia'), (30, 'Eastern Asia'), (34, 'Southern Asia'), (35, 'South-Eastern Asia'), (145, 'Western Asia'), (151, 'Eastern Europe'), (154, 'Northern Europe'), (39, 'Southern Europe'), (155, 'Western Europe'), (53, 'Australia and New Zealand'), (54, 'Melanesia'), (57, 'Micronesia'), (61, 'Polynesia')], null=True),
        ),
        migrations.AlterField(
            model_name='countryhealth',
            name='est_mdr',
            field=models.FloatField(blank=True, help_text='Estimated % Drug Resistance for the country', null=True),
        ),
    ]
