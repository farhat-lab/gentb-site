# Generated by Django 2.2.13 on 2020-12-22 13:59

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('mutations', '0023_auto_20200522_0015'),
        ('predict', '0014_auto_20190820_1130'),
    ]

    operations = [
        migrations.CreateModel(
            name='PredictResult',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('updated', models.DateTimeField(auto_now=True)),
                ('false_negative', models.DecimalField(blank=True, decimal_places=5, max_digits=10, null=True, verbose_name='False Negative Rate')),
                ('false_positive', models.DecimalField(blank=True, decimal_places=5, max_digits=10, null=True, verbose_name='False Postive Rate')),
                ('probability', models.DecimalField(blank=True, decimal_places=5, max_digits=10, null=True, verbose_name='Drug Resistance Probability')),
                ('drug', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='strain_predictions', to='mutations.Drug')),
                ('strain', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='results', to='predict.PredictStrain')),
            ],
            options={
                'unique_together': {('strain', 'drug')},
            },
        ),
        migrations.CreateModel(
            name='PredictResultLocus',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('category', models.IntegerField(choices=[(0, 'Unknown'), (1, 'Important'), (2, 'Other'), (3, 'New'), (4, 'Lineage SNPs')])),
                ('mutations', models.TextField(default='')),
                ('locus', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='prediction_results', to='mutations.GeneLocus')),
                ('result', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='loci', to='predict.PredictResult')),
            ],
            options={
                'unique_together': {('result', 'locus', 'category')},
            },
        ),
    ]
