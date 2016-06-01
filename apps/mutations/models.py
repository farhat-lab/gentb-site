import os

from django.db.models import *

from django.core.urlresolvers import reverse


class Drug(Model):
    name = CharField(max_length=255, db_index=True, unique=True)
    code = CharField(max_length=12, db_index=True, unique=True)

    class Meta:
        ordering = ('code',)

    def __str__(self):
        return '[%s] - %s' % (self.code, self.name)


class GeneLocus(Model):
    drug = ForeignKey(Drug, related_name='gene_locuses')
    name = CharField(max_length=255, db_index=True)

    class Meta:
        ordering = ('name',)
        unique_together = ('drug', 'name')

    def __str__(self):
        return self.name


class Mutation(Model):
    gene_locus = ForeignKey(GeneLocus, related_name='mutations')

    name = CharField(max_length=255, db_index=True)
    order = IntegerField(default=0)

    class Meta:
        ordering = ('order',)
        unique_together = ('gene_locus', 'name')

    def __str__(self):
        return self.name

