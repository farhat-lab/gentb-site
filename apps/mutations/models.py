#
# Copyright (C) 2016   Dr. Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# pylint: disable=too-few-public-methods, no-init, old-style-class
#
"""
Drug resistance and strain source for gene mutations django app.
"""

import os

from django.conf import settings
from django.db.models import Model, Manager, Q, QuerySet, \
    CharField, PositiveIntegerField, ForeignKey, ManyToManyField, URLField, \
    SlugField, IntegerField, BooleanField, DateField, DateTimeField, \
    TextField, DecimalField, CASCADE, SET_NULL
from django.urls import reverse

from apps.maps.models import Country, Place
from apps.uploads.models import UploadFile
from .validators import is_octal
from .utils import match_snp_name, match_snp_half, match_snp_name_raw

class DrugClassManager(Manager):
    """Allow exporting of drug classes with natural keys"""
    def get_by_natural_key(self, code):
        """The drug class code is the natural key"""
        return self.get(code=code)

class DrugClass(Model):
    """The class of drug which acts as a category for drugs"""
    name = CharField(max_length=64, db_index=True, unique=True)
    code = CharField(max_length=12, db_index=True, unique=True)

    objects = DrugClassManager()

    def natural_key(self):
        """Return the code as the lookup key for this table"""
        return (self.code,)

    def __str__(self):
        return self.name

class DrugManager(Manager):
    """Allow drugs to be exported using natural keys"""
    def get_by_natural_key(self, name):
        """The name, or the code or the abbrivation are all good keys"""
        return self.get(Q(name=name) | Q(code=name) | Q(abbr=name))

class DrugRegimen(Model):
    """A grouping of drugs that are given together in waves"""
    code = CharField(max_length=4, primary_key=True)
    name = CharField(max_length=32)
    desc = CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return self.name

class Drug(Model):
    """Each antibiotic drug which resistance might be known"""
    kind = ForeignKey(DrugClass, verbose_name='Drug Class', null=True, blank=True,
                      on_delete=SET_NULL)

    name = CharField(max_length=255, db_index=True, unique=True)
    code = CharField(max_length=12, db_index=True, unique=True)
    abbr = CharField(max_length=8, null=True, blank=True)

    priority = IntegerField(default=0, help_text="Priority of drug in regimen")
    regimen = ForeignKey(DrugRegimen, null=True, blank=True, related_name='drugs',
                         on_delete=SET_NULL)

    loci = ManyToManyField("GeneLocus", blank=True, related_name='drugs',\
        help_text="Implicated gene loci which are important to drug resistance")
    mutations = ManyToManyField("Mutation", blank=True, related_name='drugs',\
        help_text="Implicated gene mutations which cause resistance to this drug")

    objects = DrugManager()

    class Meta:
        ordering = ('regimen', '-priority',)

    def natural_key(self):
        """We want the drug code to be the key into this table"""
        return (self.code,)

    def __str__(self):
        return '[%s] - %s' % (self.code, self.name)

class GenomeManager(Manager):
    def get_by_natural_key(self, name):
        return self.get(Q(name=name) | Q(code=name))

# H37rv is a reference TBgenome, SNPs and mutations are usually called relative to this reference
class Genome(Model):
    code = SlugField(max_length=32, unique=True)
    name = CharField(max_length=255)
    length = DecimalField(default=0, max_digits=20, decimal_places=0)

    objects = GenomeManager()

    def natural_key(self):
        """The genome'sshort code is a useful key to lookup this table"""
        return (self.code,)

    def __str__(self):
        return "[%s] %s" % (self.code, self.name)


class GeneLocusManager(Manager):
    def get_by_natural_key(self, genome, name):
        return self.get(genome__code=genome, name=name)

    def for_mutation_name(self, name, brute=False):
        """Match a mutation name to a gene locus"""
        try:
            raw = match_snp_name(name)
        except ValueError:
            try:
                raw = match_snp_half(name)
            except ValueError:
                if not brute:
                    raise
                name = name.lower()
                for gene in self.all():
                    for col in ('gene_symbol', 'name', 'previous_id'):
                        if getattr(gene, col) and getattr(gene, col).lower() in name:
                            return gene
        return self._for_mutation(int(raw['ntpos']), name)

    def for_mutation(self, obj):
        """Match to a mutation object"""
        if obj.nucleotide_position is None:
            return self.for_mutation_name(obj.name)
        return self._for_mutation(obj.nucleotide_position, obj.name)

    def _for_mutation(self, pos, name):
        loci = self.filter(start__lte=pos, stop__gte=pos)
        if loci.count() == 1:
            return loci.get()
        if loci.count() == 0:
            return None

        # Overlapping loci, try and detect what we have
        (first, second) = loci
        symbol = name.replace('inter-', 'intergenic ')\
                     .replace('promoter_', 'promoter ')\
                     .split('_')[-1].lower()
        is_a = str(first.gene_symbol).lower() == symbol or str(first.name).lower() == symbol
        is_b = str(second.gene_symbol).lower() == symbol or str(second.name).lower() == symbol
        lk_a = first.name.lower() in symbol or str(first.gene_symbol).lower() in symbol
        lk_b = second.name.lower() in symbol or str(second.gene_symbol).lower() in symbol
        if (is_a and not is_b) or (lk_a and not lk_b):
            return first
        if (is_b and not is_a) or (lk_b and not lk_a):
            return second

        return None


class GeneLocus(Model):
    """
    Basically a gene as defined by a standard reference, it's location,
    length and other details which are optional.
    """
    STRANDS = (
        (None, 'Undefined'),
        ('+', '+'),
        ('-', '-'),
        ('.', '.'),
    )
    GENE_TYPES = (
        ('P', 'Promoter'),
        ('C', 'Coding'),
        ('I', 'Intergenic'),
        ('R', 'RNA'),
    )

    genome = ForeignKey(Genome, related_name='gene_locuses', null=True, blank=True,
                        on_delete=SET_NULL)
    name = CharField(max_length=255, db_index=True)
    previous_id = CharField(max_length=64, db_index=True, null=True, blank=True)

    start = IntegerField(blank=True, null=True)
    stop = IntegerField(blank=True, null=True)
    length = IntegerField(blank=True, null=True)
    strand = CharField(max_length=5, choices=STRANDS, blank=True, null=True)

    gene_type = CharField(max_length=1, choices=GENE_TYPES,\
        blank=True, null=True, help_text="Basic coding type for the locus.")
    gene_symbol = CharField(max_length=32, blank=True, null=True,\
        help_text="Short identifier used in names")
    description = CharField(max_length=255, blank=True, null=True,\
        help_text="Basic description about the gene.")

    gene_ontology = CharField(max_length=255, blank=True, null=True,\
        help_text="Gene ontology or GO-Terms are annotations in the GO format"
                  " that describe the gene product in a predictable way.")
    enzyme_commission = CharField(max_length=255, blank=True, null=True,\
        help_text="The Enzyme Commission numbers for this gene.")
    pathway_kegg = CharField(max_length=255, blank=True, null=True,\
        help_text="The KEGG based pathways")
    pathway_cyc = CharField(max_length=255, blank=True, null=True,\
        help_text="The PWY numbers usually linking to MetaCyc")
    pathway_cog = CharField(max_length=255, blank=True, null=True,\
        help_text="Clusters of Orthologous Groups of protein list")
    protein_families = CharField(max_length=255, blank=True, null=True,\
        help_text="Protein families from PFAM")

    objects = GeneLocusManager()

    class Meta:
        ordering = ('start',)
        unique_together = ('genome', 'name')

    def natural_key(self):
        """The genome's short code is a useful key to lookup this table"""
        return (self.genome.code, self.name)

    def __str__(self):
        if self.gene_symbol:
            return self.gene_symbol
        return self.name


class GeneDrugInteraction(Model):
    """
    Classify a gene drug interaction (experimental).
    """
    drug = ForeignKey(Drug, related_name='gene_interactions', on_delete=CASCADE)
    gene = ForeignKey(GeneLocus, related_name='drug_interactions', on_delete=CASCADE)
    paper = ForeignKey("Paper", related_name='interactions', null=True, blank=True,\
        help_text="Reference the paper this interaction was found in, only one blank "
                  "interaction allowed per drug/gene", on_delete=SET_NULL)
    weight = IntegerField(default=1,\
        help_text="How important is this interaction considered")
    interaction = CharField(max_length=5, default='RES', choices=[
        ('RES', 'Increases Drug Resistance'),
    ], help_text="What kind of interaction is involved in this relationship.")

    class Meta:
        unique_together = ('drug', 'gene', 'paper')

    def __str__(self):
        return "Gene '{0.gene}' {0.interaction} to '{0.drug}'".format(self)


class MutationQuerySet(QuerySet):
    @staticmethod
    def variant_names():
        """Returns a list of mutations involved in predictions"""
        # Disabled while the database and RandomForest are incongruent
        #names = self.values_list('name', flat=True)
        with open(os.path.join(settings.DATA_ROOT, 'variant_name_list.csv'), 'r') as fhl:
            return [line.strip() for line in fhl.read().split(',')]

    def matrix_csv(self, name, mutations):
        """Creates a matrix.csv file for prediction"""
        headers = ['strain']
        row = [name.replace(',', '-')]

        # Replaced with static file for now XXX:
        names = self.variant_names()

        # Mutation names could still be wrong though.
        def normalise(txt):
            for char in '!@#$%^&*_-+.,=/\\[]{}()<>\'\":;':
                txt = txt.replace(char, '_')
            return txt

        mutations = set([normalise(m) for m in mutations])

        for name in names:
            headers.append(name)
            if normalise(name) in mutations:
                row.append('1')
                mutations.remove(name)
            else:
                row.append('0')

        # we could use the csv module here, but this is ok.
        return (",".join(headers) + "\n" + ",".join(row), mutations)


class MutationManager(Manager.from_queryset(MutationQuerySet)):
    def get_by_natural_key(self, genome, locus, name):
        return self.get(gene_locus__genome__code=genome, gene_locus__name=locus, name=name)

MODE_CHOICES = (
    ("SNP", "Single Nucleotide Polymorphism"),
    ("LSP", "Long String Polymorphism"),
    ("INS", "Insertion"),
    ("DEL", "Deletion"),
)

# db_table: gtbdr.var_h37rv
class Mutation(Model):
    # genesymbol
    gene_locus = ForeignKey(GeneLocus, related_name='mutations', on_delete=CASCADE)

    # gene_id/snpname, snpid, index
    name = CharField(max_length=255, db_index=True)
    old_id = CharField(max_length=50, db_index=True, null=True, blank=True)
    order = IntegerField(default=0)

    # SNP/INS/DEL
    mode = CharField(default="SNP", max_length=4, db_index=True,\
        null=True, blank=True, choices=MODE_CHOICES)

    # ntpos, ntref, ntvar
    nucleotide_position = IntegerField(null=True, blank=True)
    nucleotide_reference = CharField(max_length=7, null=True, blank=True)
    nucleotide_varient = CharField(max_length=7, null=True, blank=True)

    # aapos, aaref, aavar
    aminoacid_position = IntegerField(null=True, blank=True)
    aminoacid_reference = CharField(max_length=41, null=True, blank=True)
    aminoacid_varient = CharField(max_length=41, null=True, blank=True,
            help_text="The variant aminoacid tracks multiple mutations "
            "in the same codon")

    # codonpos, varcodon, refcodon
    codon_position = IntegerField(null=True, blank=True)
    codon_varient = CharField(max_length=3, null=True, blank=True)
    codon_reference = CharField(max_length=3, null=True, blank=True)

    # coding, syn : [[ SNP_CN_ C and N parts ]]
    coding = CharField(max_length=1, null=True, blank=True)
    syn = CharField(max_length=1, null=True, blank=True)

    objects = MutationManager()

    def natural_key(self):
        """The genome's short code is a useful key to lookup this table"""
        return self.gene_locus.natural_key() + (self.name,)

    class Meta:
        ordering = ('order',)
        unique_together = ('gene_locus', 'name')

    def __str__(self):
        return self.name

    def save(self, **kwargs): # pylint: disable=arguments-differ
        if not self.pk and self.name:
            try:
                self.name_to_data()
            except (ValueError, KeyError):
                pass # We failed to parse the mutation
        return super().save(**kwargs)

    def name_to_data(self):
        """Unpack the name and set all the fields, return True if any of the fields changed"""
        changed = []

        _, match = match_snp_name_raw(self.name)
        data = match.groupdict()

        if 'cpos' in data and data['cpos']:
            if '-' in data['cpos']:
                data['cpos'] = data['cpos'].split('-')[-1]
            try:
                data['codpos'] = int(data['cpos']) % 3
            except ValueError:
                data['cpos'] = None

        for m_field, d_field, default in [
                ('mode', 'mode', 'SNP'),
                ('nucleotide_position', 'ntpos', None),
                ('nucleotide_reference', 'cref', None),
                ('nucleotide_varient', 'cver', None),
                ('aminoacid_position', 'apos', None),
                ('aminoacid_reference', 'aver', None),
                ('codon_position', 'codpos', None),
            ]:
            change = self._update_mutation_field(m_field, data.get(d_field, default))
            if change:
                changed.append(change)
        return changed

    def _update_mutation_field(self, m_field, data):
        """Update a sepcific field"""
        field = self._meta.get_field(m_field)
        if data:
            original = getattr(self, m_field)
            if original and not isinstance(data, type(original)):
                data = type(original)(data)

            if hasattr(field, 'max_length'):
                if field.max_length and len(data) > field.max_length:
                    raise KeyError(f"  ! {m_field} is too big len('{data}') > {field.max_length}")

            if original != data:
                setattr(self, m_field, data)
                return (m_field, original, data)
        return False



SEXES = (
    (None, 'Unknown'),
    ('F', 'Female'),
    ('M', 'Male'),
    ('O', 'Other'),
    ('', 'Left Blank'),
)
HIVS = (
    (None, 'Unknown'),
    ('negative', 'Negative'),
    ('positive', 'Positive'),
)
RESISTANCE = (
    (None, 'Unknown'),
    ('s', 'Sensitive to Drug'),
    ('i', 'Intermediate'),
    ('r', 'Resistant to Drug'),
)
RESISTANCE_GROUP = (
    (None, 'Unknown'),
    ('S', 'Sensitive'),
    ('ODR', 'Other Drug Resistant'),
    ('MDR', 'Multi Drug Resistant'),
    ('XDR', 'Extensively Drug Resistant'),
    #('TDR', 'Total Drug Resistant'),
)


# gtbdr - Targeted sequences, must link to regions table
# wgsmtb - Full sequences, must specify Null targeting.
class TargetSet(Model):
    """For targeted genes (where the whole set of mutations is imported)"""
    genome = ForeignKey(Genome, on_delete=CASCADE)
    name = CharField(max_length=64)

    def __str__(self):
        return self.name


class TargetRegion(Model):
    target_set = ForeignKey(TargetSet, related_name='regions', on_delete=CASCADE)
    gene = ForeignKey(GeneLocus, blank=True, null=True, on_delete=CASCADE)

    length = IntegerField(blank=True, null=True)
    start = IntegerField(blank=True, null=True)
    stop = IntegerField(blank=True, null=True)

    def __str__(self):
        return "Target: %s (%d-%d)" % (str(self.gene), self.start, self.stop)


class ImportSourceManager(Manager):
    """Manage collections of import sources"""
    def get_by_natural_key(self, name):
        """Get any import source by its natural key"""
        return self.get(name=name)

class ImportSource(Model):
    """Track data by how it was imported."""
    name = CharField(max_length=256)

    uploader = ForeignKey(settings.AUTH_USER_MODEL, null=True, blank=True, on_delete=SET_NULL)
    uploaded = ManyToManyField(UploadFile, blank=True, through='ImportStrain')
    complete = BooleanField(default=True)

    created = DateTimeField(auto_now_add=True)
    updated = DateTimeField(auto_now=True)

    objects = ImportSourceManager()

    def __str__(self):
        return self.name

    def natural_key(self):
        """The genome's short code is a useful key to lookup this table"""
        return (self.name,)

    def get_absolute_url(self):
        return reverse("genes:upload.view", kwargs={'pk': self.pk})

    def mutations(self):
        return StrainMutation.objects.filter(strain__importer=self)

    def resistances(self):
        if self.complete:
            return StrainResistance.objects.filter(strain__importer=self)
        return self.uploaded.get(name='resistances')

    def sources(self):
        if self.complete:
            return self.strainsource_set.all()
        return self.uploaded.get(name='sources')

    def vcf_files(self):
        """Returns a list of uploaded vcf files"""
        return self.uploaded.filter(filename__contains='vcf')

TEST_CHOICES = (
    ('', 'Not tested'),
    ('NUL', 'Value does not exist but is not required'),
    ('REQ', 'Value does not exist and is required'),
    ('NEG', 'Value is not valid'),
    ('OK', 'Value is valid'),
)

class ImportStrain(Model):
    """A link to an uploaded vcf file which we will or will have imported"""
    upload_file = ForeignKey(UploadFile, on_delete=CASCADE)
    import_source = ForeignKey(ImportSource, on_delete=CASCADE)

    is_enriched = CharField(max_length=3, choices=TEST_CHOICES, blank=True, default='')
    has_name = CharField(max_length=3, choices=TEST_CHOICES, blank=True, default='')
    has_location = CharField(max_length=3, choices=TEST_CHOICES, blank=True, default='')
    has_study = CharField(max_length=3, choices=TEST_CHOICES, blank=True, default='')
    has_patient = CharField(max_length=3, choices=TEST_CHOICES, blank=True, default='')
    has_phenotype = CharField(max_length=3, choices=TEST_CHOICES, blank=True, default='')

    def __str__(self):
        return str(self.upload_file)

class Paper(Model):
    name = CharField(max_length=128)
    doi = CharField(max_length=255, unique=True)
    url = URLField()

    notes = TextField(null=True, blank=True)

    def __str__(self):
        return self.name

class BioProject(Model):
    name = CharField(max_length=128, unique=True)

    def __str__(self):
        return self.name

class Lineage(Model):
    """
    The Coll et. al. Lineage call tree structure
    """
    slug = SlugField('Lineage ID', primary_key=True)
    name = CharField('Lineage Name', max_length=128, null=True, blank=True)
    parent = ForeignKey('self', null=True, blank=True, on_delete=CASCADE)

    def __str__(self):
        return self.name or self.slug

class StrainSourceManager(Manager):
    """Manage collections of strain sources"""
    def get_by_natural_key(self, name):
        """Get any strain source by its natural key"""
        return self.get(name=name)

class StrainSource(Model):
    """
    A strain source contains all the meta-data relating to where a strain was collected
    from whom, by whom, and what ever else we can say about the strain.
    """
    name = CharField(max_length=128, db_index=True, unique=True)
    old_id = CharField(max_length=50, db_index=True, null=True, blank=True)
    cluster = CharField(max_length=15, null=True, blank=True)
    date = DateField(null=True, blank=True)

    country = ForeignKey(Country, null=True, blank=True, related_name='sources', on_delete=SET_NULL)
    city = ForeignKey(Place, null=True, blank=True, related_name='sources', on_delete=SET_NULL)

    importer = ForeignKey(ImportSource, verbose_name='Import Source', null=True, blank=True,
                          on_delete=CASCADE)
    source_lab = CharField(max_length=100, verbose_name='Laboratory Source', db_index=True,
                           null=True, blank=True)
    source_paper = ForeignKey(Paper, related_name="strains", null=True, blank=True,
                              on_delete=SET_NULL)
    bioproject = ForeignKey(BioProject, related_name="strains", null=True, blank=True,
                            on_delete=SET_NULL)

    patient_id = CharField(max_length=16, db_index=True)
    patient_age = PositiveIntegerField(null=True, blank=True)
    patient_sex = CharField(max_length=3, choices=SEXES, null=True, blank=True)
    patient_hiv = CharField(max_length=10, choices=HIVS, null=True, blank=True)

    spoligotype_type = IntegerField(null=True, blank=True)
    spoligotype_family = CharField("Spoligotype Family Parent Strain", max_length=255, null=True, blank=True)
    spoligotype_octal = CharField(validators=[is_octal], max_length=15, null=True, blank=True)

    lineage = ForeignKey(Lineage, related_name='strains', null=True, blank=True, on_delete=SET_NULL)
    rflp_type = CharField("Restriction fragment length polymorphism type", max_length=10, null=True, blank=True)
    rflp_family = CharField("Restriction fragment length polymorphism family", max_length=10, null=True, blank=True)
    insert_type = IntegerField("Insertion sequence 6110 type", null=True, blank=True)

    wgs_group = CharField("Whole Genome Sequence Group", max_length=10, null=True, blank=True)
    principle_group = IntegerField("Principle Generic Group", null=True, blank=True)
    resistance_group = CharField(max_length=4, choices=RESISTANCE_GROUP, null=True, blank=True)
    targeting = ForeignKey(TargetSet, null=True, blank=True, on_delete=CASCADE)

    notes = TextField(null=True, blank=True)

    objects = StrainSourceManager()

    def natural_key(self):
        """Natural key for strains based just on their name"""
        return (self.name,)

    def generate_resistance_group(self):
        """Generates the correct resistance_group based on drug information"""
        resistant_to = self.drugs.filter(resistance='r').values_list('drug__code', flat=True)
        sensitive_to = self.drugs.filter(resistance='s').values_list('drug__code', flat=True)
        group = None # Unknown

        if 'INH' in resistant_to and 'RIF' in resistant_to:
            group = 'MDR'
            if ('MOXI' in resistant_to) or ('GATI' in resistant_to) or ('LEVO' in resistant_to):
                if ('KAN' in resistant_to) or ('AMK' in resistant_to) or ('CAP' in resistant_to):
                    group = 'XDR'
        elif resistant_to:
            group = 'ODR'
        elif 'INH' in sensitive_to or 'RIF' in sensitive_to:
            group = 's'

        self.resistance_group = group
        self.save()

    def __str__(self):
        return self.name or self.patient_id or ("Unnamed %d" % self.pk)

class StrainMutationManager(Manager):
    """Manage collections of mutation-strain relationships"""
    def get_by_natural_key(self, strain, genome, locus, mutation):
        """Get any mutation by natural key"""
        return self.get(strain__name=strain,
                        mutation__gene_locus__genome__code=genome,
                        mutation__gene_locus__name=locus,
                        mutation__name=mutation)

class StrainMutation(Model):
    # id (link to StrainSource imported data)
    strain = ForeignKey(StrainSource, related_name='mutations', on_delete=CASCADE)
    # snpid (via parent_id)
    mutation = ForeignKey(Mutation, related_name='strain_mutations', on_delete=CASCADE)

    # qual, depth[del], bidir
    quality = IntegerField(null=True, blank=True)

    # hqr, hqrref, fq, aavar
    mutation_reads = IntegerField(null=True, blank=True,\
        help_text="The number of sequencing reads that call the mutation")
    reference_reads = IntegerField(null=True, blank=True,\
        help_text="The number of sequencing reads that call the reference sequence")
    mapping_quality = IntegerField(null=True, blank=True,\
        help_text="Mapping quality as the root mean square of mapping qualities")
    aminoacid_varient = CharField(max_length=41, null=True, blank=True,\
        help_text="Takes into account the effect of multiple mutations in the same codon")

    objects = StrainMutationManager()

    class Meta:
        unique_together = ('strain', 'mutation')

    def natural_key(self):
        """Natural key for strains based just on their name"""
        return self.strain.natural_key() + self.mutation.natural_key()

    def __str__(self):
        return "Mutation %s for strain %s" % (str(self.mutation), str(self.strain))

    @property
    def depth(self):
        """Depth is always the mutation_reads plus reference_reads"""
        if self.mutation_reads is not None and self.reference_reads is not None:
            return self.mutation_reads + self.reference_reads

class StrainResistanceManager(Manager):
    def get_by_natural_key(self, strain, drug):
        """Get any strain resistance relationship by natural key"""
        return self.get(strain__name=strain, drug__code=drug)

class StrainResistance(Model):
    strain     = ForeignKey(StrainSource, related_name='drugs', on_delete=CASCADE)
    drug       = ForeignKey(Drug, related_name='strains', on_delete=CASCADE)
    resistance = CharField(max_length=1, choices=RESISTANCE, null=True, blank=True)

    objects = StrainResistanceManager()

    def natural_key(self):
        """Natural key for strains based just on their name"""
        return self.strain.natural_key() + self.drug.natural_key()

    def __str__(self):
        return "%s is %s to %s" % (str(self.strain), str(self.get_resistance_display()), str(self.drug))
