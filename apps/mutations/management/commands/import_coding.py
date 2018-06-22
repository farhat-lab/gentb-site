"""
Impor coding and non-coding tab files into the database
"""
import logging

from django.db import transaction, connection
from django.db.models import Q
from django.core.management.base import BaseCommand

from fav.gene_lookup import GeneLookup
from apps.mutations.models import Genome, GeneLocus

LOGGER = logging.getLogger('apps.mutations')

#connection.force_debug_cursor = True  # Change to use_debug_cursor in django < 1.8
#l = logging.getLogger('django.db.backends')
#l.setLevel(logging.DEBUG)
#l.addHandler(logging.StreamHandler())

class Command(BaseCommand):
    """Import coding command"""
    def add_arguments(self, parser):
        parser.add_argument('--files', nargs='+', type=str,\
            help='Text files to load into the mutations app')

    def handle(self, files, **options):
        """Add all the files to the coding"""
        genome = Genome.objects.get()
        gene_lookup = GeneLookup.from_files(*files)
        count = 0
        for name in gene_lookup:
            count += 1
            self.import_item(genome, name, gene_lookup[name], gene_lookup)
            if count % 100 == 0:
                print "Processed {}".format(count)
        print "finished: {}".format(count)

    @staticmethod
    def import_item(genome, name, item, lookup):
        """Impots a single item in the database"""
        try:
            # Intergenic
            before = item['geneBefore']
            after = item['geneAfter']
            if before == 'PE' == after:
                before = lookup[item['start'] - 2]['symbol'].replace('_PGRS', '')
                after = lookup[item['end'] + 2]['symbol'].replace('_PGRS', '')
            name = "{} {}-{}".format(item['type'].lower(), before, after)
            previous_id = name
            gene_type = item['type'][0].upper()
        except KeyError:
            # Genetic
            previous_id = item['synonym']
            gene_type = 'C'

        try:
            obj = GeneLocus.objects.get(name=name)
        except GeneLocus.DoesNotExist:
            try:
                obj = GeneLocus.objects.get(start=item['start'], stop=item['end'])
            except GeneLocus.DoesNotExist:
                obj = GeneLocus()

        obj.genome = genome

        obj.name = name
        obj.previous_id = previous_id
        obj.gene_type = gene_type

        obj.start = int(item['start'])
        obj.stop = int(item['end'])
        obj.strand = item['strand']
        obj.description = item['description']

        if gene_type == 'C':
            # Genetic
            obj.length = int(item['length'])
            obj.gene_symbol = item['symbol']

            obj.gene_ontology = item['goterms']
            obj.enzyme_commission = item['ecnums']
            obj.pathway_kegg = item['kegg']
            obj.pathway_cyc = item['pathways']
            obj.pathway_cog = item['cogs']
            obj.protein_families = item['pfams']
        else:
            obj.length = obj.stop - obj.start
        
        try:
            with transaction.atomic():
                obj.save()
        except Exception:
            print "Transaction ignored...."

