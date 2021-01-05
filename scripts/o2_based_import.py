#!/usr/bin/env python3
# pylint: disable=wrong-import-position

import os
import sys
import json



sys.path.insert(0, '.')
sys.path.insert(0, '..')

try:
    import manage # pylint: disable=unused-import
except ImportError as err:
    sys.stderr.write("Could not run script! Is manage.py not in the current"\
        "working directory, or is the environment not configured?:\n"\
        "{:s}\n".format(err))
    sys.exit(1)


import logging
from collections import defaultdict

from argparse import ArgumentParser
from subprocess import Popen, PIPE

from django.db import transaction
from django.db.utils import DataError

from apps.maps.models import Country, Place
from apps.maps.utils import COUNTRY_MAP, CITY_MAP

from apps.mutations.csv_lookups import Lookup as CsvLookup
from apps.mutations.models import (
    ImportSource, Genome, Lineage, Paper, Drug,
    StrainSource, BioProject, GeneLocus,
)
from apps.mutations.utils import long_match, unpack_mutation_format

LOGGER = logging.getLogger('apps.mutations')
EMPTY = {None: None, 'None': None, '': None,}
PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

class NotEnrichedError(ValueError):
    """The VCF File isn't enriched, so can't be imported"""


class Command():
    """Import VCF and JSON based strain data"""
    # Caches
    countries = EMPTY.copy()
    places = EMPTY.copy()
    studies = EMPTY.copy()
    drugs = EMPTY.copy()

    def __init__(self):
        self.arg_parser = ArgumentParser(description=self.__doc__)
        self.add_arguments(self.arg_parser)
        self.args = self.arg_parser.parse_args(sys.argv[1:])
        self.log = defaultdict(lambda: defaultdict(int))

        self.genome = Genome.objects.get(code='H37Rv')
        self.importer, created = ImportSource.objects.get_or_create(
            name=self.args.name,
            defaults={'complete': True,},
        )
        if created:
            sys.stderr.write(f"\n == Created new importer: {self.args.name} ==\n\n")
        else:
            sys.stderr.write(f"\n == Updating existing importer: {self.args.name} ==\n\n")

    def move_json(self, json_path, reason='unknown', act='reject'):
        """Reject this json, move it so it's out of the way"""
        path = os.path.join(self.args.json, f'.{act}.{reason}')
        if not os.path.isdir(path):
            os.makedirs(path)
        new_name = os.path.join(path, os.path.basename(json_path))
        os.rename(json_path, new_name)
        self.log[act][reason] += 1

    def add_arguments(self, pars):
        """Add extra command arguments"""
        pars.add_argument("-n", "--name", type=str, required=True,\
            help='The name of the import to create (or add to if an existing one)')
        pars.add_argument("-j", "--json", type=str, required=True,\
            help='Location of all the json files to import')
        pars.add_argument("-v", "--var", type=str, required=True,\
            help='The paired var files that match the json file\'s names')
        pars.add_argument("--vcf", type=str, required=True,\
            help='The paired vcf files to be turned into var files')
        return self

    def handle(self):
        """Handle the command being called"""
        for json_file in os.listdir(self.args.json):
            if not json_file.endswith('.json'):
                continue
            name = json_file.rsplit('.', 1)[0]
            var_path = os.path.join(self.args.var, name + '.var')

            if 'ID' in self.args.vcf:
                vcf_path = self.args.vcf.replace('ID', name)
            else:
                vcf_path = os.path.join(self.args.vcf, name + '.vcf')

            json_path = os.path.join(self.args.json, json_file)

            if os.path.isfile(var_path) or os.path.isfile(vcf_path):
                sys.stderr.write(f" > Importing {name}\n")

                with open(json_path, 'r') as fhl:
                    json_data = json.loads(fhl.read())

                try:
                    self.import_vcf(json_data, var_path, vcf_path)
                    self.move_json(json_path, 'done', 'ok')
                except DataError as err:
                    self.move_json(json_path, 'bad-data')
                except Country.DoesNotExist as err:
                    self.move_json(json_path, 'bad-country')
            else:
                self.move_json(json_path, 'no-vcf-or-var')

        if not self.log:
            print("== Nothing processed, no json files found ==")
        if self.log['reject']:
            print("== Rejected from Import ==")
            for reason, count in self.log['reject'].items():
                print(f" * {reason}: {count}")
        if self.log['ok']:
            print("== Imported or Ignored ==")
            for reason, count in self.log['ok'].items():
                print(f" * {reason}: {count}")

    @staticmethod
    def annotate_vcf(var_fhl, vcf_path):
        """
        Call the flat annotator and block until we're finished.
        """
        process = Popen(
            [
                os.path.join(PATH, 'bin', 'flatAnnotatorVAR_2.0.pl'),
                os.path.join(PATH, 'data', 'media', 'pipeline', 'files', 'h37rv.fasta'),
                os.path.join(PATH, 'data', 'media', 'pipeline', 'files', 'h37rv_genome_summary.txt'),
                os.path.join(PATH, 'data', 'media', 'pipeline', 'files', 'h37rv_noncoding_summary.txt'),
                vcf_path, '10', '0.4', 'PASS', 'AMB',
            ],
            shell=False, # Never have shell=True
            stdout=PIPE,
            stderr=PIPE, # Take all errors, just incase
        )
        (stdout, stderr) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write("--flat annotator failed--\n")
            sys.stderr.write(stderr.decode('utf8'))
            sys.stderr.write("--\n")
        else:
            var_fhl.write(stdout)

    @transaction.atomic
    def import_vcf(self, metadata, var_path, vcf_path):
        """Import the VCF file and list of mutations into database"""
        loc = metadata.get('location', {})
        con = loc.get('country', None)


        if not con:
            raise DataError('Location is missing from metadata')

        country = long_match(COUNTRY_MAP, self.countries, con, None, None,\
            'name', 'detail__name_short', 'detail__name_abbr', 'iso2', 'iso3',
            queryset=Country.objects.defer('geom'))

        if country is None:
            raise Country.DoesNotExist(f"Country not found: {con}")

        city = long_match(CITY_MAP, self.places, loc.get('CITY', None),\
            None, None, 'name', country=country,
            queryset=Place.objects.defer('geom'))

        pat = metadata.get('patient', {})
        datum = dict(
            importer=self.importer,
            country=country, city=city,
            patient_id=pat.get('ID', 'None'),
            patient_sex=pat.get('SEX', None),
            patient_age=(pat.get('AGE', None) or None),
            patient_hiv=pat.get('HIV', None),
            # TODO: Patient notes are discarded here, maybe keep them.
        )

        # Drug data must exist!
        drugs = metadata.get('drug', [])
        if not drugs:
            raise DataError("No phenotype drug resistance data.")

        # BIO SAMPLE first.
        name = metadata.get('sample_id', None)
        name = EMPTY.get(name, name)

        other_name = metadata.get('strain_name', None)
        other_name = EMPTY.get(other_name, other_name)

        if name is None:
            if other_name is None:
                raise DataError("No valid STRAIN_NAME or bio SAMPLE_ID")
            name, other_name = other_name, None
        datum['old_id'] = other_name

        if 'lineage' in metadata and metadata['lineage']:
            datum['lineage'] = self.get_lineage(metadata['lineage'])

        study = None
        meta_study = metadata.get('study', {})
        if 'study' in metadata:
            name = meta_study.get('name', None)
            url = meta_study.get('url', None)
            doi = meta_study.get('doi', None)
            if doi:
                study = Paper.objects.get_or_create(doi=doi, defaults={'name': name, 'url': url})[0]
            elif name:
                study = Paper.objects.get_or_create(name=name, defaults={'url': url})[0]
            if study is not None:
                datum['source_paper'] = study

        if 'project_id' in metadata:
            datum['bioproject'] = BioProject.objects.get_or_create(
                name=metadata['project_id'])[0]

        if not os.path.isfile(var_path) and os.path.isfile(vcf_path):
            sys.stderr.write(f" + Generating var {vcf_path} > {var_path}\n")
            with open(var_path, 'wb') as fhl:
                self.annotate_vcf(fhl, vcf_path)
        var = CsvLookup(filename=var_path, key='varname')

        strain, _ = StrainSource.objects.update_or_create(name=name, defaults=datum)
        for _drug in metadata.get('drug', []):
            drug_name = _drug['name']
            try:
                drug = long_match({
                    'OFLOXACIN': 'OFLX',
                }, self.drugs, drug_name, Drug, 'NOP', 'name', 'code', _match='icontains')
            except Drug.DoesNotExist:
                raise DataError(f"Can't import drug: {drug_name}")
            drug.strains.update_or_create(strain=strain, defaults=dict(
                resistance=_drug['status'].lower(),
                # TODO: More drug testing information here.
            ))

        strain.generate_resistance_group()
        self.save_mutations(strain, var)
        return name

    @staticmethod
    def save_mutations(strain, var):
        """Save all the mutations in the var file"""
        for snp in var.values():
            #gene = snp['regionid1']
            if 'varname' not in snp or len(snp['varname']) > 80:
                continue

            mutation = snp['varname']

            if len(mutation) > 150:
                sys.stderr.write("Mutation name is too large, can not add to database.\n")
                continue

            try:
                locus = GeneLocus.objects.for_mutation_name(mutation)
            except ValueError:
                continue # flatannotator outputs some whacky snp naes sometimes, ignore them.
            if locus is None:
                sys.stderr.write("Failed to unpack {varname}\n".format(**snp))
                continue

            # Correct any issues with the mutation name
            try:
                mutation = unpack_mutation_format(mutation)[2]
            except ValueError:
                continue

            try:
                (mutation, _) = locus.mutations.get_or_create(name=mutation, defaults=dict(
                    nucleotide_position=None,
                    nucleotide_reference=None,
                    nucleotide_varient=None,
                    aminoacid_position=None,
                    aminoacid_reference=None,
                    aminoacid_varient=None,
                    codon_position=snp.get('codpos', None),
                    codon_varient=snp.get('altcodon', None),
                    codon_reference=snp.get('codon', None),
                ))

                strain.mutations.get_or_create(mutation=mutation, defaults=dict(
                    quality=to_num(snp.get('qual', 0)),
                    mutation_reads=to_num(snp.get('hqr', 0)),
                    reference_reads=to_num(snp.get('hqr_ref', 0)),
                    mapping_quality=to_num(snp.get('fq', 0))))
            except (ValueError, IndexError, KeyError) as err:
                sys.stderr.write("Failed to add mutation: {} ({}) {}\n".format(mutation, err, snp))
                continue

    def get_lineage(self, name):
        parent = None
        if '.' in name:
            parent = self.get_lineage(name.rsplit('.', 1)[0])
        lineage, created = Lineage.objects.get_or_create(
            slug=name,
            defaults={'name': name, 'parent': parent})
        if not created and lineage.parent is None and parent is not None:
            lineage.parent = parent
            lineage.save()
        return lineage

def to_num(num):
    """Force a value to be a number"""
    try:
        return int(num)
    except ValueError:
        return 0

if __name__ == '__main__':
    Command().handle()
