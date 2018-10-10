
import json
from operator import or_
from collections import defaultdict

from django.core.management.base import BaseCommand, CommandError

from apps.maps.models import *
from apps.mutations.models import *
from apps.mutations.utils import *

import logging
LOGGER = logging.getLogger('apps.mutations')
EMPTY = {None: None, 'None': None, '': None,}

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--path', type=str,
            help='Directory of the data files in json format.')

    countries = EMPTY.copy()
    places = EMPTY.copy()
    MAP = {
        'unknown': None,
        'South Afica': 'South Africa',
        'Ivory Coast': 'CIV',
        'Netherland': 'Netherlands',
        'Camerun': 'Cameroon',
        'Trinidad & Tobago': 'Trinidad and Tobago',
        'Kazakstan': 'Kazakhstan',
        'Kazachstan': 'Kazakhstan',
        'Brasil': 'Brazil',
        'Azerbaidjan': 'Azerbaijan',
        'Marocco': 'Morocco',
        'DR Congo': 'Democratic Republic of the Congo',
        'RD Congo': 'Democratic Republic of the Congo',
        'Guinea Eq.': 'Equatorial Guinea',
        'Myanmar': 'Burma',
        'Netherlands Antilles': 'Sint Maarten',
        'Guinea-Conakry': 'Guinea',
        'China /Tibet': 'Tibet',
        'Carribean': 'Aruba',

        'Philipines': 'Philippines',
        'Comoro Islands': 'Comoros',
        'South Korea N': 'Korea, Republic of',
        'Korea': 'Korea, Republic of',

        'Karakalpakstan': 'Nukus', # Region
        'New York City': 'New York',
        'KwaZulu-Natal': 'Durban',
        'Kwazulu Natal': 'Durban',
        'South Carolina': 'Columbia',
        'CDC': None,
        'isolated in SF Lineage3A': None,
    }
    # Maps ids/names to Labs for the source field.
    LAB_MAP = {
        r'CDC-\d+': 'CDC',
        r'\d\d-R\d+': 'MSLI',
        r'K\d{2}\d?': 'RIVM',
        r'MT00\d\d': 'BCCDC',
        r'K-\d':         'FZB',
        r'\d{4,5}\d{2}': 'FZB',
        r'X122':        'Stellenbosch',
        r'R1207':       'Stellenbosch',
        r'[MR]\d{1,3}(Hi|L)Q': 'Stellenbosch',
        r'w-148':   'PHRI/NCBI',
        r'c':       'PHRI/NCBI',
        r'haarlem': 'PHRI/NCBI',
        r'HN878':  'NCBI',
        r'H37Ra':  'NCBI',
        r'MTB210': 'NCBI',
        r'borowski\d': 'UCSF',
        r'T\d{2}':     'UCSF',
        r'GM_0981':    'UCSF',
        r'\d{2}\d{4}': 'UCSF',
        r'M4100':      'UCSF',
    }

    def handle(self, path, **options):
        if path is None or not os.path.isdir(path):
            raise IOError("Data path doesn't exist: %s" % path)

        not_found = set(list(self.load_places(os.path.join(path, 'sources.json'), status='Loading Maps')))
        if not_found:
            sys.stderr.write('Countries or cities missing from mapping or database: %s' % str(not_found))
            return "2"

        self.make_all_drugs()

        name = os.path.basename(path.rstrip('/'))
        self.genome = Genome.objects.get(code='H37Rv')
        (self.importer, _) = ImportSource.objects.get_or_create(name=name)
        locuses = dict(self.load_genes(os.path.join(path, 'genes.json'), status='Loading Genes'))

        tarset = None
        target_file = os.path.join(path, 'targets.json')
        if os.path.isfile(target_file):
            print "Creating Genetic target set: %s" % name
            (tarset, _) = TargetSet.objects.get_or_create(name=name, genome=self.genome)
            targets = list(self.load_gene_targets(target_file, locuses=locuses, targetting=tarset, status='Loading Targets'))

        strainres = json_dict(os.path.join(path, 'resistances.json'), key_id='id')

        # We join together all the same id/snpid rows and pick the best quality one.
        mutations = json_dictlist(os.path.join(path, 'mutations.json'), key='%(id)d-%(snpid)s')
        # Quality is shown in depth.
        mutations.flatten(max, key=lambda o: o['depth'])
        # Once merged/flattened, the keys can be reset to just the id
        mutations.re_key('id')

        return list(self.load_strains(os.path.join(path, 'sources.json'), strainres, mutations, targeting=tarset, status='Loading Strains'))

    def long_match(self, d, value, model=None, *cols, **filter):
        """Match in a model with case-insensitive multi-column matching."""
        value = self.MAP.get(value, value)
        if value not in d and model and cols:
            query = reduce(or_, [Q(**{col+'__iexact': value}) for col in cols])
            d[value] = model.objects.filter(**filter).get(query)
        return d[value]

    @json_generator
    def load_places(self, row):
        """Loop through all lines and look for cities and counties"""
        try:
            country = self.long_match(self.countries, row['country'], Country, 'name', 'detail__name_short', 'detail__name_abbr', 'iso2', 'iso3')
        except Country.DoesNotExist:
            return row['country']
        try:
            self.long_match(self.places, row['city'], Place, 'name', country=country)
        except Place.DoesNotExist:
            return (row['city'], row['country'])

    @json_generator
    def load_strains(self, row, resistances, mutations, targeting=None):
        pk = str(row.pop('id'))
        res = resistances[int(pk)]
        pop_all(res, *list(row))
        _ = tr(row, ptage='patient_age', patientid='patient_id', ptsex='patient_sex', hivstatus='patient_hiv',
          snp_set=None, snp_cluster_group=None, snp_cluster_group2=None,
          sptype=('spoligotype_type', 'int'), spfamily_parentstrain='spoligotype_family', spoligo_octal='spoligotype_octal',
          inttype='rflp_type', rflpfamily='rflp_family', is6110=('insert_type', 'int'), pgg='principle_group',
          otherid=None, clustername='cluster', date=('date', 'date'), setting='source_lab', source='notes')

        row['notes'] = row.get('notes', '') or ''
        row['source_lab'] = row.get('source_lab', None)
        if row['source_lab'] and ' ' in row['source_lab']:
            row['notes'] += '; ' + row['source_lab']
            row['source_lab'] = None

        if not row['source_lab']:
            try:
                row['source_lab'] = re_match_dict(self.LAB_MAP, row['name'])
            except ValueError:
                row['source_lab'] = 'unknown'

        row['country'] = self.long_match(self.countries, row['country'])
        row['city'] = self.long_match(self.places, row['city'])
        row['resistance_group'] = res.pop('drtype')
        row['targeting'] = targeting
        row['importer'] = self.importer
        if row['name'] is None:
            row['name'] = pk
        try:
            (obj, created) = StrainSource.objects.update_or_create(defaults=row, old_id=pk)
        except Exception as err:
            print "FOUND ERR: %s" % str(err)
            return None

        for (key, v) in res.items():
            if v not in [None, '']:
                key = key[1:].upper()
                defaults = {'resistance': v}
                drug = Drug.objects.get(code=key).pk
                (dr, created) = obj.drugs.update_or_create(defaults=defaults, drug_id=drug)

        for data in mutations.get(int(pk), []):
            mutations = Mutation.objects.filter(old_id=data.pop('snpid'))
            mutation = mutations[0]
            if 'depth' in data and 'hqr' not in 'data':
                data['hqr'] = data['depth']
            _ = tr(data, qual='quality', cnsqual='quality', bidir='bi_directional', hqr='mutation_reads', hqrref='reference_reads',
              fq='mapping_quality', aavar='animoacid_varient', depth=None)
            data['bi_directional'] = data.get('bi_directional', '') == 'Y'
            (mu, created) = obj.mutations.update_or_create(defaults=data, mutation=mutation)


    @json_generator
    def load_genes(self, gene):
        gene['snpname'] = gene['snpname'].replace('.', '')
        if gene['coding'] == 'I' and gene['syn'] == 'F':
            gene['snpname'] = gene['snpname'].replace('_F_', '_CF_')
        if gene['coding'] == 'D' and gene['syn'] == 'F':
            gene['snpname'] = gene['snpname'].replace('_F_', '_CF_')
            gene['snpname'] = gene['snpname'].replace('_I_', '_CI_')
        try:
            (_, locus, mutation) = unpack_mutation_format(gene['snpname'])
        except ValueError as err:
            print "Failed to unpack %s" % str(gene)
            raise
        (locus, _) = GeneLocus.objects.update_or_create(name=locus, genome=self.genome)
        if not locus.start:
            locus.start = gene['ntpos']
            locus.previous_id = gene['geneid']
            locus.save()
        tr(gene, snpid='old_id', snpname=None, genesymbol=None, type=None, strand=None,
          description=None, geneid=None,
          ntpos='nucleotide_position', ntref='nucleotide_reference', ntvar='nucleotide_varient',
          aapos='aminoacid_position', aaref='aminoacid_reference', aavar='aminoacid_varient',
          codonpos='codon_position', varcodon='codon_varient', refcodon='codon_reference')
        (m, _) = locus.mutations.update_or_create(name=mutation, defaults=gene)
        return (gene['nucleotide_position'], locus)


    @json_generator
    def load_gene_targets(self, target, locuses, targetting=None):
        locus = None
        for location, locus in locuses.items():
            if target['stop'] > location > target['start']:
                break
        # XXX THIS WON"T WORK WITH A BREAK, SEE locus is None down there.
        if locus is None:
            if target['genesymbol']:
                gene = target['genesymbol'].replace('_', ' ')
                (locus, _) = GeneLocus.objects.get_or_create(name=gene, genome=self.genome)
                target['genesymbol'] = None
            else:
                print "Can't add %s" % str(target)
                return
        locus.gene_symbol = target.pop('genesymbol')
        locus.gene_type = target.pop('type')[0].upper()
        locus.strand = target.pop('strand')
        locus.description = target.pop('description')
        locus.save()
        target.pop('h37rv_id')
        return targetting.regions.get_or_create(gene=locus, defaults=target)


    def make_all_drugs(selselff):
        DRUGS = ['rcys', 'rmoxi', 'rpas', 'rinh', 'roflx', 'rstr', 'rkan', 'reth',
                 'ramoxclav', 'rrif', 'rlevo', 'rclof', 'rgati', 'rcap', 'rtha',
                 'rcip', 'rclar', 'rpro', 'ramk', 'rrfb', 'remb', 'rpza', 'rlin']

        drug_names = dict(
            CYS="Cycloserine",
            RFB="Rifabutin",
            CLOF="Clofazimine",
            THA="Thiacetazone",
            MOXI="Moxifloxacin",
            AMOXCLAV="Amoxicillin Clavulanate",
            GATI="Gatifloxaci",
            CLAR="Clarithromycin",
            PRO="Prothionamide",
            LIN="Linezolid",
        )

        for key in DRUGS:
            key = key[1:].upper()
            try:
                Drug.objects.get_or_create(code=key, defaults={'name': drug_names.get(key, key)})
            except Drug.DoesNotExist:
                raise ValueError("Drug '%s' doesn't exist in the database." % key)


