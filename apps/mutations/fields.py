
from django.forms.fields import CharField
from django.conf import settings


class GeneticInputField(CharField):
    def __init__(self, *args, **kwargs):
        pass

    def to_python(self, value):
        "Returns a Unicode object."
        if value in self.empty_values:
            return ''
        return smart_text(value)

    def widget_attrs(self, widget):
        attrs = super(CharField, self).widget_attrs(widget)
        if self.max_length is not None:
            # The HTML attribute is maxlength, not max_length.
            attrs.update({'maxlength': str(self.max_length)})
        return attrs

