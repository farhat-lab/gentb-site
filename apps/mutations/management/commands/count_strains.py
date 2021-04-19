"""
A command script for counting strains and populating the cache.

If you start running this command, be sure to finish or the cache will be incomplete.
"""

import sys
import logging

from collections import defaultdict

from django.db import transaction
from django.core.management.base import BaseCommand, CommandError

from apps.mutations.models import StrainMutationCache, StrainMutation, Mutation

class Command(BaseCommand):
    __doc__
    def __init__(self):
        self.rejected = 0
        self.rejects = defaultdict(int)

    def handle(self, **_):
        """Handle the command being called"""
        # Clear caching count database
        count, details = StrainMutationCache.objects.all().delete()
        print(f" [x] Deleted {count} existing caching rows.")

        total_cache = 0
        total_mut = Mutation.objects.count()
        filters = list(StrainMutationCache.matrix_filter())
        filter_names = list(StrainMutationCache.filters)
        print(f" [x] Ready to count {total_mut} Mutations.")

        for mut_x, mutation in enumerate(Mutation.objects.order_by('pk')):
            # Count all the StrainMutations
            if mut_x % 100 == 0:
                pc = mut_x / total_mut * 100
                sys.stdout.write(f" [/] Mutation: {mut_x} ({pc:0.2g}%) {total_cache} counts, {self.rejected} rejected\r")
                sys.stdout.flush()

            total_cache += self._count(mutation, filters, filter_names)

        sys.stdout.write(f" [x] Mutation: {mut_x} (100%) {total_cache} counts, {self.rejected} rejected\r")
        sys.stdout.flush()
        print("\n")

    @transaction.atomic
    def _count(self, mutation, filters, filter_names):
        """Count up all the strains linked to this mutation"""
        def _getattr(obj, name, enabled):
            if enabled:
                return None
            ret = getattr(obj, name)
            if ret is None:
                return False
            return ret

        counts = defaultdict(int)

        for stmut_x, st_mut in enumerate(mutation.strain_mutations.order_by('pk')):
            # for each combination of the given fields
            for combo in filters:
                # add one to this CacheCount
                key = tuple([_getattr(st_mut.strain, name, enabled) for name, enabled in combo])
                # Reject keys where the required field wasn't even set. This combinatin is invalid
                if False in key:
                    self.rejected += 1
                    for x, k in enumerate(key):
                        if k is False:
                            self.rejects[x] += 1
                else:
                    counts[key] += 1

        #print(f" [x] MutationStrains counted: {x} (100%) generated {ct} counts, {rejected} rejected    ")
        rejected = 0
        total = len(counts)
        #print(f" [x] Saving {ct} counts...")
        for x, (fields, count) in enumerate(counts.items()):
            #if x % 100 == 0:
                #c = x - rejected
                #pc = x / total * 100
                #sys.stdout.write(f" [ ] Caches: {c} ({pc:0.2g}%) {c} Saved {rejected} counts too small.\r")
                #sys.stdout.flush()
            if count <= 1:
                # We're not going to count single items, this is sort sorting only.
                rejected += 1
                continue

            kwargs = dict(zip(filter_names, fields))
            StrainMutationCache.objects.create(
                mutation=mutation,
                count=count,
                **kwargs)

        return total - rejected
