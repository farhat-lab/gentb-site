from django.forms import ModelForm
from apps.shared_data.models import SharedFileInfo

class SharedFileInfoForm(ModelForm):
    class Meta:
        model = SharedFileInfo
        exclude = ['md5', 'created', 'modified', 'prediction_results', 'has_prediction']
        
    def __init__(self, *args, **kwargs):
        super(SharedFileInfoForm, self).__init__(*args, **kwargs)
        # hack for bootstrap attrs
        #self.fields['contact_email'].widget.attrs.update({'type' : 'xemail'})               
        for f in self.fields.keys():
            self.fields['description'].widget.attrs.update({ 'rows' :"3"})
        
        
        #<input type="email" class="form-control" id="{{ field.id_for_label }}" placeholder="{{ field.help_text }}">
        
        