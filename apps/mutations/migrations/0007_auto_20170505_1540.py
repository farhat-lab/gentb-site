# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('mutations', '0006_strainsource_wgs_group'),
    ]

    operations = [
        migrations.AddField(
            model_name='genelocus',
            name='enzyme_commission',
            field=models.CharField(help_text=b'The Enzyme Commission numbers for this gene.', max_length=255, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='gene_ontology',
            field=models.CharField(help_text=b'Gene ontology or GO-Terms are annotations in the GO format that describe the gene product in a predictable way.', max_length=255, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='pathway_cog',
            field=models.CharField(help_text=b'Clusters of Orthologous Groups of protein list', max_length=255, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='pathway_cyc',
            field=models.CharField(help_text=b'The PWY numbers usually linking to MetaCyc', max_length=255, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='pathway_kegg',
            field=models.CharField(help_text=b'The KEGG based pathways', max_length=255, null=True, blank=True),
        ),
        migrations.AddField(
            model_name='genelocus',
            name='protein_families',
            field=models.CharField(help_text=b'Protein families from PFAM', max_length=255, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='genelocus',
            name='description',
            field=models.CharField(help_text=b'Basic description about the gene.', max_length=255, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='genelocus',
            name='gene_symbol',
            field=models.CharField(help_text=b'Short identifier used in names', max_length=32, null=True, blank=True),
        ),
        migrations.AlterField(
            model_name='genelocus',
            name='gene_type',
            field=models.CharField(blank=True, max_length=1, null=True, help_text=b'Basic coding type for the locus.', choices=[(b'P', b'Promoter'), (b'C', b'Coding'), (b'I', b'Intergenic'), (b'R', b'RNA')]),
        ),
        migrations.AlterField(
            model_name='strainsource',
            name='resistance_group',
            field=models.CharField(blank=True, max_length=4, null=True, choices=[(None, b'Unknown'), (b'S', b'Sensitive'), (b'MDR', b'Multi Drug Resistant'), (b'XDR', b'Extensively Drug Resistant')]),
        ),
    ]
