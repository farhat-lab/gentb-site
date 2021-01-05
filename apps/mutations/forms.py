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

import os

from django.conf import settings
from django.db.models import Q
from django.contrib.admin.widgets import FilteredSelectMultiple
from django.forms import (
    Form, ModelForm, CharField, Textarea, ModelMultipleChoiceField, ValidationError
)

from apps.uploads.fields import UploadField

from .models import Drug, GeneLocus, Mutation, ImportSource, ImportStrain
from .utils import unpack_mutation_format

class DrugForm(ModelForm):
    """Select drugs and mutations from a drop down"""
    paste = CharField(label="Create Mutations", widget=Textarea, required=False,\
        help_text="Put a list of mutations here for automatic processing.")
    mutations = ModelMultipleChoiceField(
        queryset=Mutation.objects.all(), label="Existing Mutations",
        widget=FilteredSelectMultiple('mutations', False))

    class Meta:
        model = Drug
        fields = ('name', 'code', 'kind', 'mutations', 'paste')

    def clean_paste(self):
        """Clean the pasted data for processing"""
        try:
            return list(self._process_data())
        except Exception as err:
            raise ValidationError('Error processing pasted data: %s' % str(err))

    def _process_data(self):
        data = self.cleaned_data.get('paste')
        for num, line in enumerate(data.strip().split("\n")):
            line = line.strip()
            if not line:
                continue
            try:
                yield unpack_mutation_format(line)
            except ValueError as err:
                raise ValueError("Line '%d:%s' can not be read: %s" % (num, line, str(err)))

    def save(self, *args, **kw):
        """Save the drug form data generating genesand mutations"""
        obj = super(DrugForm, self).save(*args, **kw)
        if obj and obj.pk:
            muts = self.cleaned_data.get('mutations')
            pks = list(muts.values_list('pk', flat=True))
            for (index, locus_name, mutation) in self.cleaned_data.get('paste'):
                locus = GeneLocus.objects.get(Q(name=locus_name) | Q(gene_symbol=locus_name))
                (mutation, _) = Mutation.objects.get_or_create(
                    name=mutation, defaults={'gene_locus': locus, 'order': index or -1})
                obj.mutations.add(mutation)
                pks.append(mutation.pk)
            self.cleaned_data['mutations'] = Mutation.objects.filter(pk__in=pks)
        return obj


class DataUploaderForm(Form):
    """
    When we want to upload new data into the maps system. This form is used.
    """
    name = CharField()
    vcf_files = UploadField(extensions=['.vcf', '.vcf.gz'], required=True, label="Enriched VCF",\
        help_text="Variant Call Formated sequence file enriched with resistance and location data.")

    @staticmethod
    def get_uploader_dir(pkey):
        """Returns the directory used to save uploaded data"""
        return os.path.join(settings.MEDIA_ROOT, 'm_uploads', str(pkey).zfill(8))

    def save(self, user):
        """
        Save the data into the importer
        """
        source = ImportSource.objects.create(
            name=self.cleaned_data.get('name'),
            uploader=user,
            complete=False,
        )
        source.save()

        path = self.get_uploader_dir(source.pk)
        if not os.path.isdir(path):
            os.makedirs(path)

        files = list(list(self.cleaned_data['vcf_files'].values())[0])
        for upload_file in files:
            upload_file.conclude_upload(path, user)
            ImportStrain.objects.create(import_source=source, upload_file=upload_file)

        source.save()
        return source
