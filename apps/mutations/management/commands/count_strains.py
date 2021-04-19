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

        def _getrkey(value):
            if value is None:
                return None
            return value is not False

        # Clear caching count database
        count, details = StrainMutationCache.objects.all().delete()
        print(f" [x] Deleted {count} existing caching rows.")

        total = StrainMutation.objects.count()
        filters = list(StrainMutationCache.matrix_filter())
        filter_names = list(StrainMutationCache.filters)
        print(f" [x] Ready to count {total} StrainMutations.")

        counts = defaultdict(int)
        rejected = 0
        rejects = defaultdict(int)

        # Count all the StrainMutations
        for y in range(0, total, 100):
            pc = (y) / total * 100
            sys.stdout.write(f" [/] StrainMutation: {y} ({pc:0.2g}%)   \r")
            sys.stdout.flush()
            for x, st_mut in enumerate(StrainMutation.objects.order_by('pk')[y:y+100]):
                if x % 50 == 5:
                    pc = (y + x) / total * 100
                    ct = len(counts)
                    sys.stdout.write(f" [\] StrainMutation: {y}+{x} ({pc:0.2g}%) {ct} counts, {rejected} rejected {rejects}      \r")
                    sys.stdout.flush()
                # for each combination of the given fields
                for combo in filters:
                    # add one to this CacheCount
                    key = tuple([_getattr(st_mut.strain, name, enabled) for name, enabled in combo])
                    # Reject keys where the required field wasn't even set. This combinatin is invalid
                    if False in key:
                        rejected += 1
                        for x, k in enumerate(key):
                            if key is False:
                                rejects[x] += 1
                    else:
                        counts[(st_mut.mutation,) + key] += 1

        print(f" [x] MutationStrains counted: {x} (100%) generated {ct} counts, {rejected} rejected    ")
        rejected = 0
        c = 0
        x = y + x + 1
        total = len(counts)
        print(f" [x] Saving {ct} counts...")
        for x, ((mutation, *fields), count) in enumerate(counts.items()):
            if x % 100 == 0:
                c = x - rejected
                pc = x / total * 100
                sys.stdout.write(f" [ ] Caches: {c} ({pc:0.2g}%) {c} Saved {rejected} counts too small.\r")
                sys.stdout.flush()
            if count <= 1:
                # We're not going to count single items, this is sort sorting only.
                rejected += 1
                continue
            kwargs = dict(zip(filter_names, fields))
            StrainMutationCache.objects.create(
                mutation=mutation,
                count=count,
                **kwargs)

        x += 1
        sys.stdout.write(f" [x] Caches: {x} (100%), {c} Saved {rejected} counts too small.\n")
