
import logging
LOGGER = logging.getLogger(__name__)

from django.views.generic import *
from django.core.urlresolvers import reverse

from .json_view import JsonView
from .models import *

class DropDownData(JsonView):
    def get_context_data(self, *kw):
        ret = {
          'levels': ['Drug', 'Gene Locus', 'Mutation'],
          'children': [],
        }
        for drug in Drug.objects.all():
            ret['children'].append({
                'name': str(drug),
                'children': [],
            })
            for locus in drug.gene_locuses.all():
                ret['children'][-1]['children'].append({
                  'name': str(locus),
                  'children': [],
                })
                for mutation in locus.mutations.all():
                    ret['children'][-1]['children'][-1]['children'].append({
                      'name': str(mutation),
                      'value': mutation.name,
                    })
        return ret

