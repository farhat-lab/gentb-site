"""
A combination of map and drug data.
"""
import json

from django.urls import reverse
from django.core.exceptions import ValidationError
from django.db.models import (
    Model, ForeignKey, SlugField, CharField, TextField, DecimalField, CASCADE, SET_NULL
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


class MapDataSource(Model):
    """A collection of rows for the data"""
    slug = SlugField(max_length=32)
    name = CharField(max_length=48)
    description = CharField(max_length=128, null=True, blank=True)

    def __str__(self):
        return self.name

    def get_columns(self):
        row = self.rows.first()
        return row is not None and row.get_columns() or []

class MapDataRow(Model):
    """A single row in the data source"""
    source = ForeignKey(MapDataSource, related_name='rows', on_delete=CASCADE)
    country = ForeignKey(Country, related_name='custom_map_rows', on_delete=CASCADE)
    drug = ForeignKey(Drug, related_name='custom_map_rows', on_delete=CASCADE)

    data = TextField(validators=[validate_json],
        help_text='The content for this row, must be a json encoded row')

    class Meta:
        ordering = ('country', 'drug')
        unique_together = (('source', 'country', 'drug'),)

    def __str__(self):
        return f"{self.parent_map}:{self.country}x{self.drug}"

    def get_data(self):
        try:
            return json.loads(self.data)
        except Exception:
            return {}

    def get_columns(self):
        return list(self.get_data().keys())

class MapDisplay(Model):
    """A displayed map of the given data"""
    slug = SlugField(max_length=32)
    name = CharField(max_length=32)
    description = CharField(max_length=128, null=True, blank=True)

    data = ForeignKey(MapDataSource, related_name='maps', on_delete=CASCADE)

    fill_column = CharField(max_length=48, default='count',
        help_text="The value used to decide what color the fill should be.")
    fill_max = DecimalField(default=-1.0, decimal_places=4, max_digits=20,
        help_text="The maximum value this field will be, if set to -1 the value will be automaticly the maximum value in the set.")
    fill_ranges = TextField(max_length=255, default='1/8,1/4,1/2,3/4',
        help_text="A comma seperated list of color ranges.")
    fill_colors = TextField(max_length=255, default='#FFFFDD,#C7E9B4,#7FCDBB,#41B6C4,#1D91C0',
        help_text="A comma seperated list of colours which the ranges relate to.")
    default_drug = ForeignKey(Drug, null=True, blank=True, on_delete=SET_NULL)

    def natural_key(self):
        return (self.slug,)

    def __str__(self):
        return self.name

    def get_max(self, values):
        fill_max = float(self.fill_max)
        if fill_max < 0.0:
            fill_max = max(values)
        return fill_max

    def get_ranges(self, fill_max):
        ret = []
        for item in self.fill_ranges.split(','):
            if '/' in item:
                a, b = [int(c) for c in item.split('/', 1)]
                item = str(float(a) / float(b))
            if '.' in item and float(item) <= 1.0:
                item = float(float(fill_max) * float(item))
            else:
                item = int(item)
            ret.append(item)
        return ret

    def get_colors(self):
        return self.fill_colors.split(',')

    def get_absolute_url(self):
        return reverse('maps:data.marginalplaces', kwargs={'slug': self.slug})

class MapDisplayDetail(Model):
    parent_map = ForeignKey(MapDisplay, related_name='details', on_delete=CASCADE)
    label = CharField(max_length=128)
    column = CharField(max_length=48, help_text="The column who's value is shown in this detail row.")
    kind = CharField(max_length=16, default='str', choices=(
        ('str', 'Text String'),
        ('float', 'Decimal Number'),
        ('int', 'Integer Number'),
    ))

    def __str__(self):
        return f"DETAILS:{self.label}"

FILTER_TYPES = (
    ('', 'No filter'),
    ('limit', 'Limit by values'),
)

class MapDisplayFilter(Model):
    parent_map = ForeignKey(MapDisplay, related_name='data_filters', on_delete=CASCADE)

    key = SlugField(max_length=32, default="")
    kind = CharField(max_length=32, choices=FILTER_TYPES, default="")
    column = CharField(max_length=48)

    label = CharField(max_length=128)
    options = TextField(help_text="A json formatted definition of the map filter options")

    def __str__(self):
        return self.label

    def get_options(self):
        try:
            return json.loads(self.options)
        except Exception as err:
            self.options = json.dumps({
                'error': str(err),
                'data': self.options,
            })
            self.save()
            return []

