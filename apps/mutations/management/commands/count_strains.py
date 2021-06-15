"""
A command script for counting strains and populating the cache.

If you start running this command, be sure to finish or the cache will be incomplete.
"""

import sys
import logging

from collections import defaultdict

from django.db import transaction
from django.db.models import Count, Q, OuterRef, Subquery
from django.core.management.base import BaseCommand, CommandError
from django.utils.timezone import now

from apps.mutations.models import (
    StrainMutationIndex, StrainMutationCount,
    StrainMutation, Mutation
)

class Command(BaseCommand):
    __doc__
    def __init__(self):
        self.rejected = 0
        self.rejects = defaultdict(int)

    def handle(self, **_):
        """Handle the command being called"""
        # Clear caching count database
        for index in StrainMutationIndex.objects.filter(generating=False, generated__isnull=True):
            index.generating = True
            index.save()
            try:
                index.count = self.generate_index(index)
            except Exception as err:
                index.delete()
                index = None
                print(f"Failed to generate_counts: '{err}'")
            if index is not None:
                index.generating = False
                index.generated = now()
                index.save()

    def generate_index(self, index):
        print("Generating index...")
        count, details = index.mutations.all().delete()
        print(f" [x] Deleted {count} existing caching rows.")

        strains = StrainMutation.objects.filter(index.strain_query('strain__')).order_by()
        count = strains.count()
        print(f" [c] Counting {count} StrainMutations into Mutations")

        strains = strains.filter(mutation=OuterRef('pk')).values('mutation')

        bulk = []
        total_count = 0
        count = 0
        for mutation_id, strain_count in Mutation.objects.annotate(
                    strain_count=Subquery(strains.annotate(c=Count('*')).values('c'))
                ).filter(
                    strain_count__gt=0,
                ).values_list('pk', 'strain_count'):
            bulk.append(
                StrainMutationCount(
                    cache_index=index,
                    mutation_id=mutation_id,
                    count=strain_count)
                )
            count += 1
            if count > 1000:
                print(f" [x] Created 1000 Counts")
                StrainMutationCount.objects.bulk_create(bulk)
                total_count += count
                count = 0
                bulk = []
        if count:
            print(f" [x] Finishing {count} Counts")
            StrainMutationCount.objects.bulk_create(bulk)
        return total_count + count
