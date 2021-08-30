"""
A combination of map and drug data.
"""
import json

from django.core.exceptions import ValidationError
from django.db.models import (
    Model, ForeignKey, SlugField, CharField, TextField, DecimalField, CASCADE,
)

from apps.maps.models import Country
from apps.mutations.models import Drug

def validate_json(text):
    """Validates submitted data as being a json row"""
    try:
        dat = json.loads(text)
        if not isinstance(dat, dict):
            raise ValidationError("Row data must be a dictionary to be valid.")
    except Exception as err:
        raise ValidationError(f"Row must be valid json: {err}")

class CustomMap(Model):
    slug = SlugField(max_length=16)
    name = CharField(max_length=32)
    description = CharField(max_length=128)

    fill_column = CharField(max_length=48, default='count',
        help_text="The value used to decide what color the fill should be.")
    fill_max = DecimalField(default=1.0, decimal_places=4, max_digits=20,
        help_text="The maximum value this field will be.")
    fill_ranges = TextField(max_length=255, default='1/8,1/4,1/2,3/4',
        help_text="A comma seperated list of color ranges.")
    fill_colors = TextField(max_length=255, default='#FFFFDD,#C7E9B4,#7FCDBB,#41B6C4,#1D91C0',
        help_text="A comma seperated list of colours which the ranges relate to.")

    def natural_key(self):
        return (self.slug,)

    def __str__(self):
        return self.name

    def get_ranges(self):
        ret = []
        for item in self.fill_ranges.split(','):
            if '/' in item:
                a, b = [int(c) for c in item.split('/', 1)]
                item = str(float(a) / float(b))
            if '.' in item and float(item) <= 1.0:
                item = float(float(self.fill_max) * float(item))
            else:
                item = int(item)
            ret.append(item)
        return ret

    def get_colors(self):
        return self.fill_colors.split(',')

class MapDetail(Model):
    parent_map = ForeignKey(CustomMap, related_name='details', on_delete=CASCADE)
    label = CharField(max_length=128)
    column = CharField(max_length=48, help_text="The column who's value is shown in this detail row.")
    kind = CharField(max_length=16, default='str', choices=(
        ('str', 'Text String'),
        ('float', 'Decimal Number'),
        ('int', 'Integer Number'),
    ))

    def __str__(self):
        return f"DETAILS:[self.label]"

class MapDataFilter(Model):
    parent_map = ForeignKey(CustomMap, related_name='data_filters', on_delete=CASCADE)

    label = CharField(max_length=128)
    column = CharField(max_length=48)
    options = TextField(help_text="A json formatted definition of the map filter options")

    def __str__(self):
        return self.label

    def get_options(self):
        try:
            return json.loads(self.options)
        except Exception:
            return []

class MapRow(Model):
    parent_map = ForeignKey(CustomMap, related_name='rows', on_delete=CASCADE)
    country = ForeignKey(Country, related_name='custom_map_rows', on_delete=CASCADE)
    drug = ForeignKey(Drug, related_name='custom_map_rows', on_delete=CASCADE)

    data = TextField(validators=[validate_json],
        help_text='The content for this row, must be a json encoded row')

    class Meta:
        ordering = ('country', 'drug')
        unique_together = (('parent_map', 'country', 'drug'),)

    def __str__(self):
        return f"{self.parent_map}:{self.country}x{self.drug}"

