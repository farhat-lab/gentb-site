# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import apps.mutations.validators


class Migration(migrations.Migration):

    dependencies = [
        ('mutations', '0001_initial'),
        ('maps', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='DrugClass',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=64, db_index=True)),
                ('code', models.CharField(unique=True, max_length=12, db_index=True)),
            ],
        ),
        migrations.CreateModel(
            name='Genome',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('code', models.SlugField(unique=True, max_length=32)),
                ('name', models.CharField(max_length=255)),
            ],
        ),
        migrations.CreateModel(
            name='StrainMutation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('quality', models.IntegerField(null=True, blank=True)),
                ('bi_directional', models.BooleanField()),
                ('mutation_reads', models.IntegerField(help_text=b'The number of sequencing reads that call the mutation', null=True, blank=True)),
                ('reference_reads', models.IntegerField(help_text=b'The number of sequencing reads that call the reference sequence', null=True, blank=True)),
                ('mapping_quality', models.IntegerField(help_text=b'Mapping quality as the root mean square of mapping qualities', null=True, blank=True)),
                ('animoacid_varient', models.CharField(help_text=b'Takes into account the effect of multiple mutations in the same codon', max_length=41, null=True, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='StrainResistance',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('resistance', models.CharField(blank=True, max_length=1, null=True, choices=[(None, b'Unknown'), (b's', b'Sensitive to Drug'), (b'i', b'Intermediate'), (b'r', b'Resistant to Drug')])),
            ],
        ),
        migrations.CreateModel(
            name='ImportSource',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=256)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('updated', models.DateTimeField(auto_now=True)),
            ],
        ),
        migrations.CreateModel(
            name='StrainSource',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=128, db_index=True)),
                ('old_id', models.CharField(db_index=True, max_length=50, null=True, blank=True)),
                ('cluster', models.CharField(max_length=15, null=True, blank=True)),
                ('date', models.DateField(null=True, blank=True)),
                ('country', models.ForeignKey(related_name='sources', blank=True, to='maps.Country', null=True)),
                ('city', models.ForeignKey(related_name='sources', blank=True, to='maps.Place', null=True)),
                ('importer', models.ForeignKey(verbose_name=b'Import Source', to='mutations.ImportSource', null=True, blank=True)),
                ('source_lab', models.CharField(max_length=100, verbose_name=b'Laboratory Source', db_index=True, null=True, blank=True)),
                ('patient_id', models.CharField(max_length=16, db_index=True)),
                ('patient_age', models.PositiveIntegerField(null=True, blank=True)),
                ('patient_sex', models.CharField(blank=True, max_length=3, null=True, choices=[(None, b'Unknown'), (b'F', b'Female'), (b'M', b'Male'), (b'O', b'Other'), (b'', b'Left Blank')])),
                ('patient_hiv', models.CharField(blank=True, max_length=10, null=True, choices=[(None, b'Unknown'), (b'negative', b'Negative'), (b'positive', b'Positive')])),
                ('spoligotype_type', models.IntegerField(null=True, blank=True)),
                ('spoligotype_family', models.CharField(max_length=255, null=True, verbose_name=b'Spoligotype Family Parent Strain', blank=True)),
                ('spoligotype_octal', models.CharField(blank=True, max_length=15, null=True, validators=[apps.mutations.validators.is_octal])),
                ('rflp_type', models.CharField(max_length=10, null=True, verbose_name=b'Restriction fragment length polymorphism type', blank=True)),
                ('rflp_family', models.CharField(max_length=10, null=True, verbose_name=b'Restriction fragment length polymorphism family', blank=True)),
                ('insert_type', models.IntegerField(null=True, verbose_name=b'Insertion sequence 6110 type', blank=True)),
                ('principle_group', models.IntegerField(null=True, verbose_name=b'Principle Generic Group', blank=True)),
                ('resistance_group', models.CharField(blank=True, max_length=4, null=True, choices=[(None, b'Unknown'), (b'S', b'Sensitive'), (b'MDR', b'Multi Drug Resistant'), (b'XDR', b'Extensively Drug Resistant'), (b'TDR', b'Total Drug Resistant')])),
                ('notes', models.TextField(null=True, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='TargetRegion',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('start', models.IntegerField(null=True, blank=True)),
                ('stop', models.IntegerField(null=True, blank=True)),
                ('length', models.IntegerField(null=True, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='TargetSet',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=64)),
                ('genome', models.ForeignKey(to='mutations.Genome')),
            ],
        ),
        migrations.AddField(
            model_name='drug',
            name='abbr',
            field=models.CharField(max_length=8, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='drug',
            name='mutations',
            field=models.ManyToManyField(help_text=b'Implicated gene mutations which cause resistance to this drug', related_name='drugs', to='mutations.Mutation', blank=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='description',
            field=models.CharField(max_length=255, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='gene_symbol',
            field=models.CharField(max_length=32, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='gene_type',
            field=models.CharField(blank=True, max_length=1, null=True, choices=[(b'P', b'Promoter'), (b'C', b'Coding'), (b'I', b'Intergenic'), (b'R', b'RNA')]),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='length',
            field=models.IntegerField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='previous_id',
            field=models.CharField(db_index=True, max_length=64, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='start',
            field=models.IntegerField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='stop',
            field=models.IntegerField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='strand',
            field=models.CharField(blank=True, max_length=5, null=True, choices=[(None, b'Undefined'), (b'+', b'+'), (b'-', b'-'), (b'.', b'.')]),
        ),
        migrations.AddField(
            model_name='mutation',
            name='aminoacid_position',
            field=models.IntegerField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='mutation',
            name='aminoacid_reference',
            field=models.CharField(max_length=41, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='mutation',
            name='aminoacid_varient',
            field=models.CharField(help_text=b'The variant aminoacid tracks multiple mutations in the same codon', max_length=41, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='mutation',
            name='coding',
            field=models.CharField(max_length=1, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='mutation',
            name='codon_position',
            field=models.IntegerField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='mutation',
            name='codon_reference',
            field=models.CharField(max_length=3, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='mutation',
            name='codon_varient',
            field=models.CharField(max_length=3, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='mutation',
            name='ecoli_aapos',
            field=models.IntegerField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='mutation',
            name='mrna_ntpos',
            field=models.IntegerField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='mutation',
            name='nucleotide_position',
            field=models.IntegerField(null=True, blank=True),
        ),
        migrations.AddField(
            model_name='mutation',
            name='nucleotide_reference',
            field=models.CharField(max_length=7, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='mutation',
            name='nucleotide_varient',
            field=models.CharField(max_length=7, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='mutation',
            name='syn',
            field=models.CharField(max_length=1, null=True, blank=True),
        ),
        migrations.AlterUniqueTogether(
            name='genelocus',
            unique_together=set([]),
        ),
        migrations.AlterUniqueTogether(
            name='mutation',
            unique_together=set([]),
        ),
        migrations.AddField(
            model_name='targetregion',
            name='gene',
            field=models.ForeignKey(blank=True, to='mutations.GeneLocus', null=True),
        ),
        migrations.AddField(
            model_name='targetregion',
            name='target_set',
            field=models.ForeignKey(related_name='regions', to='mutations.TargetSet'),
        ),
        migrations.AddField(
            model_name='strainsource',
            name='targeting',
            field=models.ForeignKey(blank=True, to='mutations.TargetSet', null=True),
        ),
        migrations.AddField(
            model_name='strainresistance',
            name='drug',
            field=models.ForeignKey(related_name='strains', to='mutations.Drug'),
        ),
        migrations.AddField(
            model_name='strainresistance',
            name='strain',
            field=models.ForeignKey(related_name='drugs', to='mutations.StrainSource'),
        ),
        migrations.AddField(
            model_name='strainmutation',
            name='mutation',
            field=models.ForeignKey(related_name='strain_mutations', to='mutations.Mutation'),
        ),
        migrations.AddField(
            model_name='strainmutation',
            name='strain',
            field=models.ForeignKey(related_name='mutations', to='mutations.StrainSource'),
        ),
        migrations.AddField(
            model_name='drug',
            name='kind',
            field=models.ForeignKey(verbose_name=b'Drug Class', blank=True, to='mutations.DrugClass', null=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='genome',
            field=models.ForeignKey(related_name='gene_locuses', blank=True, to='mutations.Genome', null=True),
        ),
    ]
