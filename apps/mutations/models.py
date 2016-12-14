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

class MutationQuerySet(QuerySet):
    def matrix_csv(self, name, mutations):
        """Creates a matrix.csv file for prediction"""
        mutations = set(mutations)
        headers = ['strain']
        row = [name.replace(',', '-')]
        for item in self:
            headers.append(item.name)
            if item.name in mutations:
                row.append('1')
                mutations.remove(item.name)
            else:
                row.append('0')

        # we could use the csv module here, but this is ok.
        return (",".join(headers) + "\n" + ",".join(row), mutations)

class Mutation(Model):
    gene_locus = ForeignKey(GeneLocus, related_name='mutations')

    name = CharField(max_length=255, db_index=True)
    order = IntegerField(default=0)
    objects = MutationQuerySet.as_manager()

    class Meta:
        ordering = ('order',)
        unique_together = ('gene_locus', 'name')

    def __str__(self):
        return self.name

