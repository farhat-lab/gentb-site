from django import forms
from apps.predict.models import PredictDataset,\
                PredictDatasetStatus,\
                DATASET_STATUS_UPLOADED_READY_ID


# requires loading of apps/predict/fixtures/initial_data.json
#


class UploadPredictionDataForm(forms.ModelForm):

    """error_messages = {
        'password_mismatch': _("The two password fields didn't match."),
        'password_complexity_length': _("The password must be at least 7 characters long (and contain at least 1 letter and 1 digit.)"),
        'password_complexity_content': _("The password must contain at least 1 letter and 1 digit."),
    }

    affiliation = forms.CharField(label='Affiliation', max_length=255)
    retype_password = forms.CharField(label='Retype Password', widget=forms.PasswordInput())
    """
    class Meta:
        model = PredictDataset
        fields = ('title', 'description', 'file1', 'file2')

    def get_vfc_dataset(self, tb_user):

        assert tb_user is not None, "tb_user cannot be None"
        assert hasattr(self, 'cleaned_data'), "Do not call this method on an invalid form."

        # save PredictDataset, made inactive
        #
        vcf_dataset = self.save(commit=False)   # get object
        vcf_dataset.user = tb_user              # set user
        vcf_dataset.has_prediction = False
        vcf_dataset.set_status_uploaded_ready(save_status=False)
        vcf_dataset.save()      # save the object

        return vcf_dataset


class DatasetRunNotificationForm(forms.Form):

    run_md5 = forms.CharField()
    success = forms.BooleanField(required=False)
    result_data = forms.CharField(widget=forms.Textarea, required=False)

    def was_run_successful(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['success']

    def get_run_md5(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['run_md5']

    def get_result_data(self):
        assert hasattr(self, 'cleaned_data'), "Must have is_valid() == true.  cleaned_data not found"
        return self.cleaned_data['result_data']


