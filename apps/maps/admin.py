from django.contrib.admin import *

from .models import Place, Country

site.register(Place)
site.register(Country)

