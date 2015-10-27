from django import forms
from apps.predict.models import PredictDataset,\
                PredictDatasetStatus
                #DATASET_STATUS_UPLOADED_READY_ID


# requires loading of apps/predict/fixtures/initial_data.json
#


class UploadPredictionDataForm(forms.ModelForm):
    """
    Form for a user to enter a title, description, and dropbox_url

    The dropbox_url is used to retrieve dropbox metadata
    """
    class Meta:
        model = PredictDataset
        fields = ('title', 'description', 'dropbox_url', )

    def get_dataset(self, tb_user):

        assert hasattr(self, 'cleaned_data'), "Do not call this method on an invalid form. (call is_valid() first)"
        assert tb_user is not None, "tb_user cannot be None"

        # -------------------------------------
        # save PredictDataset, made inactive
        # -------------------------------------
        ds = self.save(commit=False)   # get object
        ds.user = tb_user              # set user
        ds.set_status_not_ready(save_status=False)
        ds.save()      # save the object

        return ds


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
