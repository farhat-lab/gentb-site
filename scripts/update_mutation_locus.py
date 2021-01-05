#!/usr/bin/env python3
# pylint: disable=wrong-import-position

import sys

sys.path.insert(0, '.')
sys.path.insert(0, '..')

try:
    import manage # pylint: disable=unused-import
except ImportError as err:
    sys.stderr.write("Could not run script! Is manage.py not in the current"\
        "working directory, or is the environment not configured?:\n"\
        "{:s}\n".format(err))
    sys.exit(1)

from django.db.models import Q
from apps.utils.models import merge_model_objects
from apps.mutations.models import GeneLocus, Mutation
from apps.mutations.utils import match_snp_name, match_snp_half

GENES = GeneLocus.objects.filter(start__isnull=False, stop__isnull=False)

def update_mutation_info(mut, mode='SNP', syn=None, ntpos=None, cref=None, cpos=None, cver=None,
                         aref=None, apos=None, aver=None, noncode=None, codes=None, **other):
    if mode.upper() in ('SNP', 'LSP', 'INS', 'DEL'):
        mut.mode = mode.upper()

    try:
        mut.nucleotide_position = int(ntpos)
    except (ValueError, TypeError):
        pass

    if cref:
        mut.nucleotide_reference = cref[:7]

    if cver:
        mut.nucleotide_varient = cver[:7]

    if aref:
        mut.aminoacid_reference = aref[:41]
    if aver:
        mut.aminoacid_varient = aver[:41]

    try:
        mut.aminoacid_position = int(apos)
    except (ValueError, TypeError):
        pass

    try:
        mut.codon_position = int(cpos) % 3
    except (ValueError, TypeError):
        pass

    # codon_varient = 
    # codon_reference = 

    if syn and syn[0] in 'PICN':
        mut.coding = syn[0]
        if len(syn) > 1 and syn[1] in 'DIFSZ':
            mut.syn = syn[1]

def clean_mutation(mut):
    """Clean the mutation information as much as possible"""
    # Set other values
    try:
        raw = match_snp_name(mut.name)
    except ValueError:
        print(f"Can't clean: {mut.name}")
    else:
        update_mutation_info(mut, **raw)
        sys.stdout.write(f" * {mut.name}: {mut.nucleotide_position} {mut.nucleotide_reference}->{mut.nucleotide_varient} A:{mut.aminoacid_reference}->{mut.aminoacid_varient}          \r")
        mut.save()

def set_gene_locus(mut):
    others = Mutation.objects.filter(name=mut.name).exclude(pk=mut.pk)
    if others.count() > 0:
        print("Merging with an existing mutation {}".format(mut))
        merge_model_objects(mut, list(others), keep_old=True)
        total, counts = others.delete()
        counts.pop('mutations.Mutation', None)
        if set(counts.values()) != {0}:
            print("Deleted some unexpected values when merging mutations: {}".format(counts))
            sys.exit(1)

    if mut.nucleotide_position is None:
        try:
            raw = match_snp_name(mut.name)
        except ValueError:
            raw = match_snp_half(mut.name)
        mut.nucleotide_position = int(raw['ntpos'])
    
    # Set the right gene locus
    loci = GENES.filter(start__lte=mut.nucleotide_position, stop__gte=mut.nucleotide_position)
    if loci.count() == 1:
        locus = loci.get()
        if locus != mut.gene_locus:
            print("Moving mutation from {} to {}".format(mut.gene_locus, locus))
            mut.gene_locus = locus
        else:
            print("Already OK.")
    elif loci.count() == 0:
        print("No gene or intergenic for mutation {}".format(mut))
    else:
        (a, b) = loci
        symbol = str(mut).replace('inter-', 'intergenic ')\
                         .replace('promoter_', 'promoter ')\
                         .split('_')[-1].lower()
        is_a = str(a.gene_symbol).lower() == symbol or str(a.name).lower() == symbol
        is_b = str(b.gene_symbol).lower() == symbol or str(b.name).lower() == symbol
        if is_a and not is_b:
            print("A: Moving mutation from {} to {}".format(mut.gene_locus, a))
            mut.gene_locus = a
        elif is_b and not is_a:
            print("B: Moving mutation from {} to {}".format(mut.gene_locus, a))
            mut.gene_locus = b
        else:
            if a.name.lower() in symbol or str(a.gene_symbol).lower() in symbol:
                print("C: Moving mutation from {} to {}".format(mut.gene_locus, a))
                mut.gene_locus = a
            elif b.name.lower() in symbol or str(b.gene_symbol).lower() in symbol:
                print("D: Moving mutation from {} to {}".format(mut.gene_locus, b))
                mut.gene_locus = b
            else:
                print("NO: {2} {3} == {0.name} ({0.gene_symbol}) vs. {1.name} ({1.gene_symbol})".format(a, b, mut, symbol))

    mut.save()

if __name__ == '__main__':
    # Make sure every mutation name can be processed into the correct data.
    #failed = []
    #for mut in Mutation.objects.all():
    #    try:
    #        raw = match_snp_name(mut.name)
    #    except ValueError:
    #        failed.append(mut.name)
    #    if len(failed) == 40:
    #        print("Mercy rules, 40 failed mutation names:")
    #        break

    #if failed:
    #    print("\n-- Some mutation names can't be processed --\n")
    #    for item in failed:
    #        print(item)
    #    sys.exit(0)

    for mut in Mutation.objects.filter(gene_locus__stop__isnull=True):
        set_gene_locus(mut)

    for mut in Mutation.objects.filter(Q(nucleotide_position__isnull=True) | Q(mode__isnull=True)):
        clean_mutation(mut)

    gloc = GeneLocus.objects.filter(mutations__isnull=True, stop__isnull=True)
    total, counts = gloc.delete()
    if total:
        print("Deleted {} bad loci: {}".format(total, counts))

    print("Loci remaining: {}".format(GeneLocus.objects.all().count()))
