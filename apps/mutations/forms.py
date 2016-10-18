
from django.forms import ModelForm, CharField, Textarea, ValidationError

from .models import Drug
from .utils import unpack_mutation_format

class DrugForm(ModelForm):
    paste = CharField(label="Paste Mutations", widget=Textarea, required=False,
               help_text="Put a list of mutations here for automatic processing.")

    class Meta:
        model = Drug
        fields = ('name', 'code', 'paste')

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
            for (index, locus_name, mutation) in self.cleaned_data.get('paste'):
                (loc, _) = obj.gene_locuses.get_or_create(name=locus_name)
                (mut, _) = loc.mutations.get_or_create(name=mutation,
                        defaults={'order': index or -1})
        return obj
