
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

class Command(BaseCommand):
    # Caches
    countries = EMPTY.copy()
    places = EMPTY.copy()

    def handle(self, **kw):
        try:
            self.pipeline = Pipeline.objects.get(name='VCF_TO_VAR')
        except Pipeline.DoesNotExist:
            sys.stderr.write("No pipeline 'VCF_TO_VAR' which is required.\n")
            return

        for importer in ImportSource.objects.filter(complete=False, uploader__isnull=False):
            try:
                self.import_source(importer)
            except UploadFile.DoesNotExist:
                sys.stderr.write(" [!] Missing upload file for importer: %s\n" % str(importer))
            except (FileNotFound, FieldsNotFound) as err:
                sys.stderr.write(" [!] %s: %s\n" % (str(err), str(importer)))

    def import_source(self, importer):
        count = importer.vcf_files().count()
        uploads = importer.vcf_files().filter(retrieval_end__isnull=False)
        notloads = count - uploads.count()
        if notloads:
            return self.importer.set_status("Waiting for %d Uploads", notloads)

        # Create a pipeline to process all the uploads
        for upload in uploads:
            if upload.flag not in ['VCF', 'ERROR']:
                runner = self.pipeline.run(
                    'IMPORTER%d:VCF%d' % (importer.pk, upload.pk),
                    file=upload.fullpath)
                if runner.update_all():
                    upload.flag = 'VCF'
                    upload.save()

        ready = uploads.filter(flag='VCF')
        notready = count - ready.count()
        if notready:
            return importer.set_status("Waiting for %d VCF Files", notready)

        self.genome = Genome.objects.get(code='H37Rv')
        self.importer = importer

        for fl in self.importer.vcf_files():
            vcf = VCFReader(filename=fl.fullpath)
            var = CsvLookup(filename=fl.fullpath[:-4] + '.var', key='varname')
            try:
                self.import_vcf(vcf, var)
            except ValueError as err:
                fl.flag = 'ERR'
                fl.error = str(err)
                fl.save()
            except Country.DoesNotExist:
                fl.flag = 'NOCO'
                fl.save()
            

    def import_vcf(self, vcf, var):
        if 'LOCATION' not in vcf.metadata:
            raise ValueError('NO_METADATA')

        loc = vcf.metadata['LOCATION'][0]
        country = long_match(COUNTRY_MAP, self.countries, loc.get('COUNTRY', None),
            Country, 'name', 'detail__name_short', 'detail__name_abbr', 'iso2', 'iso3')

        city = long_match(CITY_MAP, self.places, loc.get('CITY', None), Place, None, 'name', country=country)

        pat = vcf.metadata.get('PATIENT', [{}])[0]
        strain = StrainSource.objects.get_or_create(name=vcf.metadata['STRAIN_NAME'], defaults=dict(
            country=country, city=city,
            patient_id=pat.get('ID', None),
            patient_sex=pat.get('SEX', None),
            patient_age=pat.get('AGE', None),
            patient_hiv=pat.get('HIV', None),
            # TODO: Patient notes are discarded here, maybe keep them.
        ))

        study = long_march({}, self.studies, vcf.metadata['STUDY'].get('NAME', None), 'name', 'url', default=None)

        for drug in vcf.metadata.get('DRUG', []):
            drug = Drug.objects.get(Q(name__icontains=drug['NAME']) | Q(code__icontains=drug['NAME']))
            drug.strains.update_or_create(strain=strain, defaults=dict(
                resistance=drug['STATUS'],
                # TODO: More drug testing information here.
            ))

        for snp in var:
            gene = snp['regionid1']
            try:
                (_, locus, mutation) = unpack_mutation_format(snp['varname'])
                # All genes in the gene summary should be already loaded.
                locus = GeneLocus.objects.get(name=locus, genome=self.genome)
            except ValueError as err:
                raise ValueError("Failed to unpack {varname}".format(**snp))
            except GeneLocus.DoesNotExist:
                raise ValueError("Failed to get gene {varname}, are all genes loaded from reference?".format(**snp))

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
                codon_varient=None,
                codon_reference=None,
            ))

