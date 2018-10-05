
import sys
import json
import logging

from vcf import VCFReader

from collections import defaultdict

from django.core.management.base import BaseCommand, CommandError

from apps.maps.models import Country, Place
from apps.maps.utils import COUNTRY_MAP, CITY_MAP
from apps.uploads.models import UploadFile
from apps.pipeline.models import Pipeline

from apps.mutations.csv_lookups import Lookup as CsvLookup
from apps.mutations.models import *
from apps.mutations.utils import *

LOGGER = logging.getLogger('apps.mutations')
EMPTY = {None: None, 'None': None, '': None,}

class DataError(ValueError):
    pass

class NotEnrichedError(ValueError):
    pass


class Command(BaseCommand):
    # Caches
    countries = EMPTY.copy()
    places = EMPTY.copy()
    studies = EMPTY.copy()
    drugs = EMPTY.copy()

    def handle(self, **kw):
        try:
            self.pipeline = Pipeline.objects.get(name='VCF_TO_VAR')
        except Pipeline.DoesNotExist:
            sys.stderr.write("No pipeline 'VCF_TO_VAR' which is required.\n")
            return

        for importer in ImportSource.objects.filter(complete=False, uploader__isnull=False):
            try:
                if not self.import_source(importer):
                    # Only break if we returned a 'WAIT' signal
                    continue
            except UploadFile.DoesNotExist:
                sys.stderr.write(" [!] Missing upload file for importer: {}\n".format(str(importer)))
            except (FileNotFound, FieldsNotFound, DataError) as err:
                sys.stderr.write(" [!] {}: {}\n".format(str(err), str(importer)))

            importer.complete = True
            importer.save()

    def import_source(self, importer):
        count = importer.vcf_files().count()
        uploads = importer.vcf_files().filter(retrieval_end__isnull=False)
        notloads = count - uploads.count()
        if notloads:
            return sys.stderr.write("Waiting for {} Uploads\n".format(notloads))

        uploads.filter(flag='ERR').update(flag='OK', retrieval_error='')

        # Create a pipeline to process all the uploads
        for upload in uploads.exclude(flag='VCF'):
            output_dir = os.path.dirname(upload.fullpath)
            runner = self.pipeline.run(
                'IMPORTER{:d}:VCF{:d}'.format(importer.pk, upload.pk),
                rerun=True, file=upload.fullpath, output_dir=output_dir)
            if upload.flag not in ['VCF', 'ERR']:
                if runner.update_all():
                    if runner.get_errors():
                        upload.flag = 'ERR'
                        upload.retrieval_error = runner.get_errors()
                    else:
                        upload.flag = 'VCF'
                    upload.save()

        err = uploads.filter(flag='ERR')
        if err.count():
            raise DataError("Errors in VAR loading caused us to stop.")

        for upload in uploads.exclude(flag='VCF'):
            if os.path.isfile(upload.fullpath[:-4] + '.var'):
                print("Manual flag {} as DONE".format(upload.fullpath))
                upload.flag = 'VCF'
                upload.save()
            
        ready = uploads.filter(flag='VCF')
        notready = count - ready.count()
        if err.count() == 0 and notready:
            print(uploads.exclude(flag='VCF'))
            return sys.stderr.write("Waiting for {} VCF Files\n".format(notready))

        self.genome = Genome.objects.get(code='H37Rv')

        for fl in importer.vcf_files():
            fl.retrieval_error = ""
            vcf = VCFReader(filename=fl.fullpath)
            var_file = fl.fullpath[:-4] + '.var'
            if not os.path.isfile(var_file):
                print("Failed to find var file: {}".format(var_file))
                continue
            var = CsvLookup(filename=fl.fullpath[:-4] + '.var', key='varname')
            try:
                fl.name = self.import_vcf(importer, vcf, var)
                fl.flag = 'DONE'
            except DataError as err:
                fl.flag = 'DATA'
                sys.stderr.write("Error: {}\n".format(err))
                fl.retrieval_error = str(err)
            except NotEnrichedError:
                fl.flag = 'POOR'
            except Country.DoesNotExist as err:
                fl.flag = 'NOCO'
                fl.retrieval_error = str(err)

            fl.save()
        return True

    def import_vcf(self, importer, vcf, var):
        """Import the VCF file and list of mutations into database"""
        if 'ENRICHED' not in vcf.metadata:
            raise NotEnrichedError()

        if 'LOCATION' not in vcf.metadata:
            raise DataError('Location is missing from metadata')

        loc = vcf.metadata['LOCATION'][0]
        con = loc.get('COUNTRY', None)
        country = long_match(COUNTRY_MAP, self.countries, con, Country, None,\
            'name', 'detail__name_short', 'detail__name_abbr', 'iso2', 'iso3')

        if country is None:
            raise DataError("Country not found: {}".format(con))

        city = long_match(CITY_MAP, self.places, loc.get('CITY', None),\
            Place, None, 'name', country=country)

        pat = vcf.metadata.get('PATIENT', [{}])[0]
        datum = dict(
            importer=importer,
            country=country, city=city,
            patient_id=pat.get('ID', None),
            patient_sex=pat.get('SEX', None),
            patient_age=(pat.get('AGE', None) or None),
            patient_hiv=pat.get('HIV', None),
            # TODO: Patient notes are discarded here, maybe keep them.
        )

        # Drug data must exist!
        drugs = vcf.metadata.get('DRUG', [])
        if not drugs:
            raise DataError("No phenotype drug resistance data.")

        # BIO SAMPLE first.
        name = vcf.metadata.get('SAMPLE_ID', (None,))[0]
        name = EMPTY.get(name, name)

        other_name = vcf.metadata.get('STRAIN_NAME', (None,))[0]
        other_name = EMPTY.get(other_name, other_name)

        if name is None:
            if other_name is None:
                raise DataError("No valid STRAIN_NAME or bio SAMPLE_ID")
            name, other_name = other_name, None

        datum['old_id'] = other_name
        strain, _ = StrainSource.objects.update_or_create(name=name, defaults=datum)

        #study = None
        #if 'STUDY' in vcf.metadata:
        #    study = long_match({}, self.studies, vcf.metadata['STUDY']\
        #        .get('NAME', None), Paper, 'name', 'url', 'doi', default=None)

        for _drug in vcf.metadata.get('DRUG', []):
            try:
                drug = long_match({
                    'OFLOXACIN': 'OFLX',
                }, self.drugs, _drug['NAME'], Drug, 'NOP', 'name', 'code', _match='icontains')
            except Drug.DoesNotExist:
                raise DataError("Can't import drug: {}".format(_drug['NAME']))
            drug.strains.update_or_create(strain=strain, defaults=dict(
                resistance=_drug['STATUS'].lower(),
                # TODO: More drug testing information here.
            ))

        strain.generate_resistance_group()

        for snp in var.values():
            #gene = snp['regionid1']
            if len(snp['varname']) > 80:
                continue
            try:
                (_, locus, mutation) = unpack_mutation_format(snp['varname'])
                # All genes in the gene summary should be already loaded.
                #locus = GeneLocus.objects.get(name=locus, genome=self.genome)
                # Ignore that and ask for a new one
                locus, _ = GeneLocus.objects.update_or_create(name=locus, genome=self.genome)
            except ValueError:
                #raise DataError("Failed to unpack {varname}".format(**snp))
                sys.stderr.write("Failed to unpack {varname}\n".format(**snp))
                continue
            except GeneLocus.DoesNotExist:
                raise DataError("Failed to get gene {varname}, "\
                    "are all genes loaded from reference?".format(**snp))

            if len(mutation) > 150:
                sys.stderr.write("Mutation name is too large, can not add to database.\n")
                continue

            locus.mutations.get_or_create(name=mutation, defaults=dict(
                nucleotide_position=None,
                nucleotide_reference=None,
                nucleotide_varient=None,
                aminoacid_position=None,
                aminoacid_reference=None,
                aminoacid_varient=None,
                codon_position=snp['codpos'],
                codon_varient=snp['altcodon'],
                codon_reference=snp['codon'],
            ))

        return name
