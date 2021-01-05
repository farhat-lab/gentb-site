# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from collections import defaultdict
from django.db import migrations

def to_drug_mutations(apps, schema_editor):
    """Move the drug<->mutation relationship to a ManyToMany field and remove duplicate mutations"""
    Mutation = apps.get_model("mutations", "Mutation")
    GeneLocus = apps.get_model("mutations", "GeneLocus")
    Drug = apps.get_model("mutations", "Drug")
    Genome = apps.get_model("mutations", "Genome")

    # All mutations so far entered have been h37rv
    (genome, _) = Genome.objects.get_or_create(name="Mycobacterium tuberculosis", code="H37Rv")

    drug_mutations = defaultdict(set)
    for pkg in Drug.objects.values('pk', 'gene_locuses__mutations__name'):
        drug_mutations[pkg['pk']].add(pkg['gene_locuses__mutations__name'])

    unique_locus = {}
    for locus in GeneLocus.objects.all():
        if locus.name in unique_locus:
            # Add all of these mutations to the master locus
            unique_locus[locus.name].mutations.add(*locus.mutations.all())
            # Now get rid of this pretender!
            locus.delete()
        else:
            # Record this unique locus using it's name as the unique property
            unique_locus[locus.name] = locus
            locus.genome = genome
            locus.save()

    unique_mutations = {}
    for mutation in Mutation.objects.all():
        if mutation.name in unique_mutations:
            # Delete this duplicate mutation name
            mutation.delete()
        else:
            # Record this unique mutation using it's name as above
            unique_mutations[mutation.name] = mutation

    # Now we need to reattach drugs to mutations
    for (drug_pk, mutations) in drug_mutations.items():
        drug = Drug.objects.get(pk=drug_pk)
        # Add a list of mutation objects based on their names
        drug.mutations.add(*[unique_mutations[name] for name in mutations])


def to_drug_locuses(apps, schema_editor):
    """Move the relationship back to a foreignkey field and create duplicate mutations"""
    Mutation = apps.get_model("mutations", "Mutation")
    GeneLocus = apps.get_model("mutations", "GeneLocus")
    Drug = apps.get_model("mutations", "Drug")

    for mutation in Mutation.objects.all():
        for x, drug in enumerate(mutation.drugs.all()):
            if x == 0:
                # For the first drug, we set the drug directly to the gene_locus
                # Basically using the same object without making a new one.
                mutation.gene_locus.drug = drug
            else:
                # For all drugs not the first one, we make a new gene_locus and new
                # mutation object and chain them from the drug. (you can tell right away
                # why this structure was a bad idea and why we moved away from it)
                new_locus = GeneLocus.objects.create(name=mutation.gene_locus.name, drug=drug)
                new_mutation = Mutation.objects.create(name=mutation.name, gene_locus=new_locus)

    GeneLocus.objects.filter(drug__isnull=True).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('mutations', '0002_auto_20161118_1041'),
    ]

    operations = [
        migrations.RunPython(to_drug_mutations, to_drug_locuses),
    ]
