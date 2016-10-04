
from django.forms import ModelForm, CharField, Textarea, ValidationError

from .models import Drug

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
                yield self._process_line(line)
            except ValueError as err:
                raise ValueError("Line '%d:%s' can not be read: %s" % (no, line, str(err)))

    def _process_line(self, line):
        index = None
        if " " in line:
            index, line = line.split(" ", 1)
            try:
                index = int(index)
            except:
                raise ValueError("Optional sort index should be a number.")

        bits = line.split('_')
        if bits[1] == 'P':
            return (index, 'promoter ' + bits[-1], line)
        elif bits[1] == 'I':
            return (index, 'intergenic ' + bits[-1], line)
        elif bits[1] in ['CN', 'CD', 'CF', 'CI', 'CZ', 'N', 'ND', 'NI', 'NF']:
            return (index, bits[-1], line)

        raise ValueError("Must be promoter, intergenic or CN, CD, CF, CI, CZ or N, ND, NI, NF")

    def save(self, *args, **kw):
        obj = super(DrugForm, self).save(*args, **kw)
        if obj and obj.pk:
            for (index, locus_name, mutation) in self.cleaned_data.get('paste'):
                (loc, _) = obj.gene_locuses.get_or_create(name=locus_name)
                (mut, _) = loc.mutations.get_or_create(name=mutation,
                        defaults={'order': index or -1})
        return obj
