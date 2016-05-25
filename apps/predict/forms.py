"""
Django forms for adding PredictDataset objects as well as a Confirmation form

"""
import json

from django.forms import * 

from apps.utils.file_patterns import FILE_TYPE_FASTQ, FilePatternHelper
from apps.dropbox.models import DropboxRetrievalLog
from apps.dropbox.util import get_dropbox_metadata_from_link
from apps.dropbox.widgets import DropboxChooserWidget
from apps.mutations.fields import GeneticInputField

from .models import PredictDataset

class UploadConfirmForm(ModelForm):
    """Confirm the files in the upload"""
    class Meta:
        model = PredictDataset
        fields = ('status',)

class ManualInputForm(ModelForm):
    """
    Manually enter genetic information for prediction.
    """
    genetic_information = GeneticInputField()

    class Meta:
        model = PredictDataset
        fields = ('title', 'description')


class UploadForm(ModelForm):
    """
    Form for a user to enter a title, description, and dropbox_url

    The dropbox_url is used to retrieve dropbox metadata
    """
    class Meta:
        model = PredictDataset
        exclude = ('md5', 'file_directory', 'has_prediction')
        widgets = {
            'dropbox_url': Textarea(attrs={'rows': '4'}),
            'status': HiddenInput(),
            'user': HiddenInput(),
            #'dropbox_url': DropboxChooserWidget(['.fastq', '.vcf']),
        }

    def clean_fastq_type(self):
        """
        If this is a FastQ file, make sure the user has chosen a FastQ type
        """
        file_type = self.cleaned_data.get('file_type', None)
        fastq_type = self.cleaned_data.get('fastq_type')
        if file_type == FILE_TYPE_FASTQ and not fastq_type:
            raise ValidationError("For FastQ files, please choose a FastQ type")
        return fastq_type

    def clean_dropbox_url(self):
        """Check the dropbox metadata from the url link"""
        # (This should be moved to an async or ajax call in another part of the code )
        url = self.cleaned_data.get('dropbox_url')
        file_type = self.cleaned_data.get('file_type', None)
        file_patterns = FilePatternHelper.get_file_patterns_for_dropbox(file_type)

        # Use the dropbox API to look at the files under this dropbox_url
        (success, self.box) = get_dropbox_metadata_from_link(url, file_patterns=file_patterns)

        if not success:
            raise ValidationError(self.box)

        return url

    def save(self, **kw):
        dataset = super(type(self), self).save(**kw)
        if dataset.pk:
            DropboxRetrievalLog.objects.create(dataset=dataset,
                file_metadata=self.box.dropbox_link_metadata,
                selected_files=self.box.matching_files_metadata,
            )
        return dataset


class NotificationForm(Form):
    success = BooleanField(required=False)
    result_data = CharField(widget=Textarea, required=False)

    def was_run_successful(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['success']

    def get_result_data(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['result_data']

