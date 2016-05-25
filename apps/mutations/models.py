import os

from django.db.models import *

from django.core.urlresolvers import reverse


class Drug(Model):
    name = CharField(max_length=255, db_index=True)
    code = CharField(max_length=12, db_index=True)

    class Meta:
        ordering = ('code',)

    def __str__(self):
        return '[%s] - %s' % (self.code, self.name)


class GeneLocus(Model):
    drug = ForeignKey(Drug, related_name='genes')
    name = CharField(max_length=255, db_index=True)

    class Meta:
        ordering = ('name',)

    def __str__(self):
        return self.name


class Mutation(Model):
    gene_locus = ForeignKey(GeneLocus, related_name='mutations')

    name = CharField(max_length=255, db_index=True)
    order = IntegerField(default=0)

    class Meta:
        ordering = ('order',)

    def __str__(self):
        return self.name

