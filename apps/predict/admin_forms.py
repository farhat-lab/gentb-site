from django import forms
import os
from os.path import isdir

class PipelineScriptsDirectoryForm(forms.ModelForm):

    def clean_script_directory(self):
        script_dir = self.cleaned_data['script_directory']

        script_dir = script_dir.strip()
        if len(script_dir) < 4:
            raise forms.ValidationError('Please give a valid directory name (more than 4 chars): %s' % script_dir)

        if not isdir(script_dir):
            raise forms.ValidationError('This directory does not exist: %s' % script_dir)

        return script_dir
