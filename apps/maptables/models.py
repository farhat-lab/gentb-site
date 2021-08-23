"""
A combination of map and drug data.
"""
import json

from django.core.exceptions import ValidationError
from django.db.models import (
    Model, ForeignKey, SlugField, CharField, TextField, CASCADE,
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

    def natural_key(self):
        return (self.slug,)

    def __str__(self):
        return self.name

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

