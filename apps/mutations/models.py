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
"""
Drug resistance and strain source for gene mutations django app.
"""

import os

from django.conf import settings
from django.db.models import *
from django.core.urlresolvers import reverse

from .validators import is_octal
from apps.maps.models import Country, Place

class DrugClass(Model):
    name = CharField(max_length=64, db_index=True, unique=True)
    code = CharField(max_length=12, db_index=True, unique=True)

    def __str__(self):
        return self.name


class Drug(Model):
    kind = ForeignKey(DrugClass, verbose_name='Drug Class', null=True, blank=True)

    name = CharField(max_length=255, db_index=True, unique=True)
    code = CharField(max_length=12, db_index=True, unique=True)
    abbr = CharField(max_length=8, null=True, blank=True)

    mutations = ManyToManyField("Mutation", blank=True, related_name='drugs',
        help_text="Implicated gene mutations which cause resistance to this drug")

    class Meta:
        ordering = ('code',)

    def __str__(self):
        return '[%s] - %s' % (self.code, self.name)


# H37rv is a reference TBgenome, SNPs and mutations are usually called relative to this reference
class Genome(Model):
    code = SlugField(max_length=32, unique=True)
    name = CharField(max_length=255)

    def __str__(self):
        return "[%s] %s" % (self.code, self.name)


class GeneLocus(Model):
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

    genome = ForeignKey(Genome, related_name='gene_locuses', null=True, blank=True)
    name = CharField(max_length=255, db_index=True)
    previous_id = CharField(max_length=64, db_index=True, null=True, blank=True)

    start  = IntegerField(blank=True, null=True)
    stop   = IntegerField(blank=True, null=True)
    length = IntegerField(blank=True, null=True)
    strand = CharField(max_length=5, choices=STRANDS, blank=True, null=True)

    gene_type   = CharField(max_length=1, choices=GENE_TYPES,
        blank=True, null=True, help_text="Basic coding type for the locus.")
    gene_symbol = CharField(max_length=32, blank=True, null=True,
        help_text="Short identifier used in names")
    description = CharField(max_length=255, blank=True, null=True,
        help_text="Basic description about the gene.")

    gene_ontology = CharField(max_length=255, blank=True, null=True,
        help_text="Gene ontology or GO-Terms are annotations in the GO format"
        " that describe the gene product in a predictable way.")
    enzyme_commission = CharField(max_length=255, blank=True, null=True,
        help_text="The Enzyme Commission numbers for this gene.")
    pathway_kegg = CharField(max_length=255, blank=True, null=True,
        help_text="The KEGG based pathways")
    pathway_cyc = CharField(max_length=255, blank=True, null=True,
        help_text="The PWY numbers usually linking to MetaCyc")
    pathway_cog = CharField(max_length=255, blank=True, null=True,
        help_text="Clusters of Orthologous Groups of protein list")
    protein_families = CharField(max_length=255, blank=True, null=True,
        help_text="Protein families from PFAM")

    class Meta:
        ordering = ('name',)
        unique_together = ('genome', 'name')

    def __str__(self):
        return self.name


class MutationQuerySet(QuerySet):
    def matrix_csv(self, name, mutations):
        """Creates a matrix.csv file for prediction"""
        headers = ['strain']
        row = [name.replace(',', '-')]

        # Disabled while the database and RandomForest are incongruent
        #names = self.values_list('name', flat=True)
        # Replaced with static file for now XXX:
        fn = os.path.join(settings.DATA_ROOT, 'variant_name_list.csv')
        with open(fn, 'r') as fhl:
            names = [line.strip() for line in fhl.read().split(',')]

        # Mutation names could still be wrong though.
        def normalise(txt):
            for char in '!@#$%^&*_-+.,=/\\[]{}()<>\'\":;':
                txt = txt.replace(char, '_')
            return txt

        #mutations = set(mutations)
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
    pass


# db_table: gtbdr.var_h37rv
class Mutation(Model):
    # genesymbol
    gene_locus = ForeignKey(GeneLocus, related_name='mutations')

    # gene_id/snpname, snpid, index
    name   = CharField(max_length=255, db_index=True)
    old_id = CharField(max_length=50, db_index=True, null=True, blank=True)
    order  = IntegerField(default=0)

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

    mrna_ntpos = IntegerField(null=True, blank=True)
    ecoli_aapos = IntegerField(null=True, blank=True)

    predictor = BooleanField(default=False,
        help_text="This mutation is selected to be used in predictions and "
        "will be shown to users in the manual mutation selection process.")

    objects = MutationManager()

    class Meta:
        ordering = ('order',)
        unique_together = ('gene_locus', 'name')

    def __str__(self):
        return self.name

SEXES = (
  (None, 'Unknown'),
  ('F', 'Female'),
  ('M', 'Male'),
  ('O', 'Other'),
  ('',  'Left Blank'),
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
   ('MDR', 'Multi Drug Resistant'),
   ('XDR', 'Extensively Drug Resistant'),
   #('TDR', 'Total Drug Resistant'),
)


# gtbdr - Targeted sequences, must link to regions table
# wgsmtb - Full sequences, must specify Null targeting.
class TargetSet(Model):
    genome = ForeignKey(Genome)
    name = CharField(max_length=64)

    def __str__(self):
        return self.name


class TargetRegion(Model):
    target_set = ForeignKey(TargetSet, related_name='regions')
    gene       = ForeignKey(GeneLocus, blank=True, null=True)

    start       = IntegerField(blank=True, null=True)
    stop        = IntegerField(blank=True, null=True)
    length      = IntegerField(blank=True, null=True)

    def __str__(self):
        return "Target: %s (%d-%d)" % (str(self.gene), self.start, self.stop)


class ImportSource(Model):
    """Track data by how it was imported."""
    name = CharField(max_length=256)

    created = DateTimeField(auto_now_add=True)
    updated = DateTimeField(auto_now=True)

    def __str__(self):
        return self.name


# db_table: gtbdr.Strainsourcedata
class StrainSource(Model):
    name    = CharField(max_length=128, db_index=True)
    old_id  = CharField(max_length=50, db_index=True, null=True, blank=True)
    cluster = CharField(max_length=15, null=True, blank=True)
    date    = DateField(null=True, blank=True)

    country = ForeignKey(Country, null=True, blank=True, related_name='sources')
    city    = ForeignKey(Place, null=True, blank=True, related_name='sources')

    importer = ForeignKey(ImportSource, verbose_name='Import Source', null=True, blank=True)
    source_lab = CharField(max_length=100, verbose_name='Laboratory Source', db_index=True, null=True, blank=True)

    patient_id  = CharField(max_length=16, db_index=True)
    patient_age = PositiveIntegerField(null=True, blank=True)
    patient_sex = CharField(max_length=3, choices=SEXES, null=True, blank=True)
    patient_hiv = CharField(max_length=10, choices=HIVS, null=True, blank=True)

    spoligotype_type   = IntegerField(null=True, blank=True)
    spoligotype_family = CharField("Spoligotype Family Parent Strain", max_length=255, null=True, blank=True)
    spoligotype_octal  = CharField(validators=[is_octal], max_length=15, null=True, blank=True)

    rflp_type       = CharField("Restriction fragment length polymorphism type", max_length=10, null=True, blank=True)
    rflp_family     = CharField("Restriction fragment length polymorphism family", max_length=10, null=True, blank=True)
    insert_type     = IntegerField("Insertion sequence 6110 type", null=True, blank=True)

    wgs_group = CharField("Whole Gnome Sequence Group", max_length=10, null=True, blank=True)
    principle_group = IntegerField("Principle Generic Group", null=True, blank=True)
    resistance_group = CharField(max_length=4, choices=RESISTANCE_GROUP, null=True, blank=True)
    targeting = ForeignKey(TargetSet, null=True, blank=True)

    notes = TextField(null=True, blank=True)

    def __str__(self):
        return self.name or self.patient_id or ("Unnamed %d" % self.pk)


class StrainMutation(Model):
    # id (link to StrainSource imported data)
    strain = ForeignKey(StrainSource, related_name='mutations')
    # snpid (via parent_id)
    mutation = ForeignKey(Mutation, related_name='strain_mutations')

    # qual, depth[del], bidir
    quality = IntegerField(null=True, blank=True)
    bi_directional = BooleanField()

    # hqr, hqrref, fq, aavar
    mutation_reads    = IntegerField(null=True, blank=True, help_text="The number of sequencing reads that call the mutation")
    reference_reads   = IntegerField(null=True, blank=True, help_text="The number of sequencing reads that call the reference sequence")
    mapping_quality   = IntegerField(null=True, blank=True, help_text="Mapping quality as the root mean square of mapping qualities")
    animoacid_varient = CharField(max_length=41, null=True, blank=True, help_text="Takes into account the effect of multiple mutations in the same codon")

    def __str__(self):
        return "Mutation %s for strain %s" % (str(self.mutation), str(self.strain))

    @property
    def depth(self):
        """Depth is always the mutation_reads plus reference_reads"""
        if self.mutation_reads is not None and self.reference_reads is not None:
            return self.mutation_reads + self.reference_reads


class StrainResistance(Model):
    strain     = ForeignKey(StrainSource, related_name='drugs')
    drug       = ForeignKey(Drug, related_name='strains')
    resistance = CharField(max_length=1, choices=RESISTANCE, null=True, blank=True)

    def __str__(self):
        return "%s is %s to %s" % (str(self.strain), str(self.get_resistance_display()), str(self.drug))

