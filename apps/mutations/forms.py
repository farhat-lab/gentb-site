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
Drop mutation for pasting in mutation information quickly.
"""

from django.forms import *
from django.contrib.admin.widgets import FilteredSelectMultiple

from .models import Drug, GeneLocus, Mutation
from .utils import unpack_mutation_format

class DrugForm(ModelForm):
    paste = CharField(label="Create Mutations", widget=Textarea, required=False,
        help_text="Put a list of mutations here for automatic processing.")
    mutations = ModelMultipleChoiceField(
        queryset=Mutation.objects.all(), label="Existing Mutations",
        widget=FilteredSelectMultiple('mutations', False))

    class Meta:
        model = Drug
        fields = ('name', 'code', 'mutations', 'paste')

    def clean_paste(self):
        try:
            return list(self._process_data())
        except Exception as err:
            raise ValidationError('Error processing pasted data: %s' % str(err))

    def _process_data(self):
        data = self.cleaned_data.get('paste')
        for no, line in enumerate(data.strip().split("\n")):
            line = line.strip()
            if not line:
                continue
            try:
                yield unpack_mutation_format(line)
            except ValueError as err:
                raise ValueError("Line '%d:%s' can not be read: %s" % (no, line, str(err)))

    def save(self, *args, **kw):
        obj = super(DrugForm, self).save(*args, **kw)
        if obj and obj.pk:
            muts = self.cleaned_data.get('mutations')
            pks = list(muts.values_list('pk', flat=True))
            for (index, locus_name, mutation) in self.cleaned_data.get('paste'):
                (locus, _) = GeneLocus.objects.get_or_create(name=locus_name)
                (mutation, _) = Mutation.objects.get_or_create(name=mutation,
                    defaults={'gene_locus': locus, 'order': index or -1})
                obj.mutations.add(mutation)
                pks.append(mutation.pk)
            self.cleaned_data['mutations'] = Mutation.objects.filter(pk__in=pks)
        return obj
