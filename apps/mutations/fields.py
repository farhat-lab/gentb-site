
from django.forms.fields import CharField
from django.forms.widgets import Textarea
from django.conf import settings

class ListerWidget(Textarea):
    class Media:
        js = ('js/lister.js',)

class GeneticInputField(CharField):
    def __init__(self, data_url, *args, **kw):
        self.data_url = data_url
        kw['widget'] = ListerWidget
        super(GeneticInputField, self).__init__(*args, **kw)

    def to_python(self, value):
        "Returns a Unicode object."
        if value in self.empty_values:
            return ''
        return smart_text(value)

    def widget_attrs(self, widget):
        attrs = super(CharField, self).widget_attrs(widget)
        attrs['data-data-url'] = self.data_url
        attrs['class'] = 'lister'
        return attrs

