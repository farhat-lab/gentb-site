"""
A command script for counting strains and populating the cache.

If you start running this command, be sure to finish or the cache will be incomplete.
"""

import sys
import logging

from collections import defaultdict
from django.core.management.base import BaseCommand, CommandError

from apps.mutations.models import StrainMutationCache, StrainMutation

class Command(BaseCommand):
    __doc__
    def handle(self, **_):
        """Handle the command being called"""
        def _getattr(obj, name, enabled):
            if enabled:
                return None
            ret = getattr(obj, name)
            if ret is None:
                return False
            return ret

        # Clear caching count database
        count, details = StrainMutationCache.objects.all().delete()
        print(f"Deleted {count} existing caching rows.")

        total = StrainMutation.objects.count()
        filters = StrainMutationCache.matrix_filter()
        filter_names = list(StrainMutationCache.filters)

        counts = defaultdict(int)

        # Count all the StrainMutations
        for x, st_mut in enumerate(StrainMutation.objects.all()):
            if x % 200 == 0:
                pc = x / total * 100
                ct = len(counts)
                sys.stdout.write(f"StrainMutation: {x} ({pc}%) generated {ct} counts so far\r")
            # for each combination of the given fields
            for combo in filters:
                # add one to this CacheCount
                key = tuple([_getattr(st_mut.strain, name, enabled) for name, enabled in combo])
                # Reject keys where the required field wasn't even set. This combinatin is invalid
                if False in key:
                    continue
                counts[(st_mut.mutation,) + key] += 1

        total = len(counts)
        print("\nSaving counts...\n")
        for x, ((mutation, *fields), count) in enumerate(counts.items()):
            if x % 200 == 0:
                pc = x / total * 100
                print(f"Caches Saved: {x} ({pc}%)")
            kwargs = dict(zip(filter_names, fields))
            StrainMutationCache.objects.create(
                mutation=mutation,
                count=count,
                **kwargs)

