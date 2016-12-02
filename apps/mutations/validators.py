from django.core.exceptions import ValidationError
from django.utils.translation import ugettext_lazy as _

def is_octal(value):
    """Validate that an incoming value is an octal value"""
    value = str(value)
    if not value.isdigit() or "8" in value or "9" in value:
        raise ValidationError(_('%(value)s is not an octal number'), params={'value': value})

