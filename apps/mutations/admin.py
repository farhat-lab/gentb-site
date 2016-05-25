from django.contrib.admin import *
from .models import *

site.register(Drug)
site.register(GeneLocus)
site.register(Mutation)

