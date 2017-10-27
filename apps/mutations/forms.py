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

from django.forms import *
from django.contrib.admin.widgets import FilteredSelectMultiple

from apps.uploads.models import UploadFile

from .models import Drug, GeneLocus, Mutation, ImportSource
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


from apps.uploads.fields import UploadField, UploadTable

from django.conf import settings
def get_uploader_dir(pk):
    return os.path.join(settings.MEDIA_ROOT, 'm_uploads', str(pk).zfill(8))

class DataUploaderForm(Form):
    """
    When we want to upload new data into the maps system. This form is used.
    """
    name = CharField()
    sources = UploadTable(required=True,
        help_text="Each source should contain the country and city fields so"
              " they can be placed on the map.",
        columns=['id',
        "ptage", "city", "otherid", "inttype", "name", "rflpfamily", "sptype",
        "clustername", "country", "hivstatus", "patientid",
        "spfamily_parentstrain", "source", "setting", "ptsex", "is6110",
        "date", "spoligo_octal", "pgg"],
    )
    vcf_files = UploadField(extensions=['.vcf', '.vcf.gz'], required=True,
        label="Source VCF", help_text="Variant Call Formated sequence data"\
        " file. Each file should be named with the [source_id].vcf or when"\
        " ordered should match the order in the sources csv file.")

    resistances = UploadTable(required=True,
        columns=['source_id', 'drug', 'resistant'],
        parsers=['\d+', '\w+', '[01]'])

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

        path = get_uploader_dir(source.pk)
        if not os.path.isdir(path):
            os.makedirs(path)

        files = list(self.cleaned_data['vcf_files'].values()[0])
        for upload_file in files:
            upload_file.conclude_upload(path, user)

        # Create a file for the tb_source csv and resistances list csv
        tb_sources = UploadFile.objects.create(
            name='sources', filename='sources.csv',
            file_directory=path)
        tb_sources.save_now(self.cleaned_data['sources'])
        files.append(tb_sources)

        resistances = UploadFile.objects.create(
            name='resistances', filename='resistances.csv',
            file_directory=path)
        resistances.save_now(self.cleaned_data['resistances'])
        files.append(resistances)

        source.uploaded = files
        source.save()
        return source


