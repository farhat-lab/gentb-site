
import sys
import logging

from vcf import VCFReader

from django.db import transaction
from django.db.utils import DataError
from django.core.management.base import BaseCommand, CommandError

from apps.maps.models import Country, Place
from apps.maps.utils import COUNTRY_MAP, CITY_MAP
from apps.uploads.models import UploadFile

from apps.mutations.csv_lookups import Lookup as CsvLookup
from apps.mutations.models import (
    ImportSource, Genome, Lineage, Paper, Drug,
    StrainSource, BioProject, GeneLocus,
    StrainMutation,
)
from apps.mutations.utils import *

LOGGER = logging.getLogger('apps.mutations')
EMPTY = {None: None, 'None': None, '': None,}

class NotEnrichedError(ValueError):
    """The VCF File isn't enriched, so can't be imported"""


class Command(BaseCommand):
    # Caches
    countries = EMPTY.copy()
    places = EMPTY.copy()
    studies = EMPTY.copy()
    drugs = EMPTY.copy()

    def __init__(self):
        self.genome = Genome.objects.get(code='H37Rv')

    def handle(self, **_):
        """Handle the command being called"""
        for importer in ImportSource.objects.filter(complete=False, uploader__isnull=False):
            done = 0
            try:
                for uploaded_file in importer.vcf_files().filter(retrieval_end__isnull=False):
                    done += self.import_mutations(importer, uploaded_file)
                    uploaded_file.save()
            except UploadFile.DoesNotExist:
                sys.stderr.write(" [!] Missing upload file: {}\n".format(str(importer)))
            except (FileNotFound, FieldsNotFound, DataError) as err:
                sys.stderr.write(" [!] {}: {}\n".format(str(err), str(importer)))

            importer.complete = (done == importer.vcf_files().count())
            importer.save()

    def import_mutations(self, importer, vcf_file):
        """
        Import mutations from this vcf file, returns True if the process is complete
        either in error or success. retrieval_error is set when in error while flag
        is set when waiting.
        """
        if vcf_file.flag not in ["", None, "WAIT"]:
            # Previously processed, ignore
            return True

        vcf_file.retrieval_error = ""
        if not os.path.isfile(vcf_file.fullpath):
            # File is missing, even though uploader says the uploading process is complete.
            vcf_file.retrieval_error = "VCF File is missing!"
            vcf_file.flag = 'MIA'
            return True

        var_file = vcf_file.fullpath[:-4] + '.var'
        if not os.path.isfile(var_file):
            # It's not ready yet, an external process will populate the var file
            vcf_file.flag = "WAIT"
            return False

        vcf = VCFReader(filename=vcf_file.fullpath)
        var = CsvLookup(filename=vcf_file.fullpath[:-4] + '.var', key='varname')
        try:
            vcf_file.name = self.import_vcf(importer, vcf, var)
            vcf_file.flag = 'DONE'
        except DataError as err:
            vcf_file.flag = 'DATA'
            sys.stderr.write("Error: {}\n".format(err))
            vcf_file.retrieval_error = str(err)
        except NotEnrichedError:
            vcf_file.flag = 'POOR'
        except Country.DoesNotExist as err:
            vcf_file.flag = 'NOCO'
            vcf_file.retrieval_error = str(err)
        return True

    @transaction.atomic
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
            patient_id=pat.get('ID', 'None'),
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

        if 'LINEAGE' in vcf.metadata:
            datum['lineage'] = Lineage.objects.get_or_create(
                name=vcf.metadata['LINEAGE'][0],
                slug=name,
            )[0]

        study = None
        if 'STUDY' in vcf.metadata:
            name = vcf.metadata['STUDY'][0].get('NAME', None)
            url = vcf.metadata['STUDY'][0].get('URL', None)
            doi = vcf.metadata['STUDY'][0].get('DOI', None)
            if doi:
                study = Paper.objects.get_or_create(doi=doi, defaults={'name': name, 'url': url})[0]
            elif name:
                study = Paper.objects.get_or_create(name=name, defaults={'url': url})[0]
            if study is not None:
                datum['source_paper'] = study

        if 'PROJECT_ID'in vcf.metadata:
            datum['bioproject'] = BioProject.objects.get_or_create(
                name=vcf.metadata['PROJECT_ID'][0])[0]

        strain, _ = StrainSource.objects.update_or_create(name=name, defaults=datum)

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
        self.save_mutations(strain, var)
        return name

    def save_mutations(self, strain, var):
        """Save all the mutations in the var file"""
        for snp in var.values():
            #gene = snp['regionid1']
            if 'varname' not in snp or len(snp['varname']) > 80:
                continue
            try:
                (_, locus, mutation) = unpack_mutation_format(snp['varname'])
                # All genes in the gene summary should be already loaded.
                locus = GeneLocus.objects.get(name=locus, genome=self.genome)
            except ValueError:
                sys.stderr.write("Failed to unpack {varname}\n".format(**snp))
                continue
            except GeneLocus.DoesNotExist:
                raise DataError("Failed to get gene {varname}, "\
                    "are all genes loaded from reference?".format(**snp))

            if len(mutation) > 150:
                sys.stderr.write("Mutation name is too large, can not add to database.\n")
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
                    mapping_quality=to_num(snp.get('fq', 0)))
                )
            except (ValueError, IndexError, KeyError) as err:
                sys.stderr.write("Failed to add mutation: {} ({}) {}\n".format(mutation, err, snp))
                continue

def to_num(num):
    """Force a value to be a number"""
    try:
        return int(num)
    except ValueError:
        return 0
