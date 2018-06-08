"""
Impor coding and non-coding tab files into the database
"""
import logging

from django.core.management.base import BaseCommand

from fav.gene_lookup import GeneLookup
from apps.mutations.models import Genome, GeneLocus

LOGGER = logging.getLogger('apps.mutations')

class Command(BaseCommand):
    """Import coding command"""
    def add_arguments(self, parser):
        parser.add_argument('--files', nargs='+', type=str,\
            help='Text files to load into the mutations app')

    def handle(self, files, **options):
        """Add all the files to the coding"""
        genome = Genome.objects.get()
        gene_lookup = GeneLookup.from_files(*files)
        for name in gene_lookup:
            self.import_item(genome, name, gene_lookup[name])

    @staticmethod
    def import_item(genome, name, item):
        """Impots a single item in the database"""
        try:
            obj = GeneLocus.objects.get(name=name)
        except GeneLocus.DoesNotExist:
            try:
                obj = GeneLocus.objects.get(start=item['start'], stop=item['end'])
            except GeneLocus.DoesNotExist:
                obj = GeneLocus()

        obj.genome = genome
        obj.name = name

        obj.start = item['start']
        obj.stop = item['end']
        obj.strand = item['strand']
        obj.description = item['description']

        try:
            # Intergenic
            obj.gene_type = item['type'][0].upper()
            obj.previous_id = item['geneBefore'] + item['geneAfter']
        except KeyError:
            # Genetic
            obj.length = item['length']
            obj.previous_id = item['synonym']
            obj.gene_symbol = item['symbol']

            obj.gene_ontology = item['goterms']
            obj.enzyme_commission = item['ecnums']
            obj.pathway_kegg = item['kegg']
            obj.pathway_cyc = item['pathways']
            obj.pathway_cog = item['cogs']
            obj.protein_families = item['pfams']
        obj.save()
