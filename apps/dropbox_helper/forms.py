"""
Form used to validate callback parameters
"""
from os.path import isdir
import json

from django import forms
from django.core.validators import RegexValidator

class DropboxRetrievalParamsForm(forms.Form):

    dropbox_url = forms.URLField()
    destination_directory = forms.CharField()
    callback_url = forms.URLField()
    callback_md5 = forms.CharField(validators=[RegexValidator(\
                    regex='^\w{32}$',\
                    message='Not a valid MD5. Must be 32 characters.',\
                    code='invalid_md5'\
                )],)

    def clean_destination_directory(self):
        dest_dir = self.cleaned_data.get('destination_directory', None)

        if not isdir(dest_dir):
            raise forms.ValidationError("Destination directory not found: {0}".format(dest_dir))

        return dest_dir

    def get_data_as_json_string(self):
        d = dict(dropbox_url=self.cleaned_data['dropbox_url'],\
                destination_directory=self.cleaned_data['destination_directory'],\
                callback_url=self.cleaned_data['callback_url'],\
                callback_md5=self.cleaned_data['callback_md5'],\
                )
        return json.dumps(d)
