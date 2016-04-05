"""Add a template tag to turn python objects into JSON"""
import types
import json

from django import template
from django.utils.safestring import mark_safe

register = template.Library()

@register.filter
def jsonify(obj):
    if isinstance(obj, types.GeneratorType):
        obj = list(obj)
    return mark_safe(json.dumps(obj))

