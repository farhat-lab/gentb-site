"""
A command script for counting strains and populating the cache.

If you start running this command, be sure to finish or the cache will be incomplete.
"""

import sys
import logging

from time import perf_counter
from contextlib import contextmanager
from collections import defaultdict

from django.db import transaction
from django.db.models import Count, Q, OuterRef, Subquery
from django.core.management.base import BaseCommand, CommandError
from django.utils.timezone import now

from apps.mutations.models import (
    StrainMutationIndex, StrainMutationCount,
    StrainSource, StrainMutation, Mutation
)

from django.core.paginator import Paginator

@contextmanager
def timer():
    start = perf_counter()
    yield lambda: perf_counter() - start


def chunked_iterator(queryset, chunk_size=10000):
    paginator = Paginator(queryset, chunk_size)
    for page in range(1, paginator.num_pages + 1):
        for obj in paginator.page(page).object_list:
            yield (page / paginator.num_pages) * 100, obj

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
            #index.save()
            try:
                index.count = self.generate_index(index)
            except Exception as err:
                #index.delete()
                index = None
                print(f"\n [!] Failed to generate_counts: '{err}'\n")
            if index is not None:
                index.generating = False
                index.generated = now()
                #index.save()

    def _generate_index(self, index):
        print("Generating index...")
        count, details = index.mutations.all().delete()
        print(f" [x] Deleted {count} existing caching rows.")

        strains = StrainSource.objects.filter(index.strain_query(''))
        count = strains.count()
        print(f" [c] Counting {count} Strains into Mutations")

        bulk = []
        bulk_count = 0
        st_count = 0
        mut_count = 0
        this_count = 0
        this_mutation = None
        threshold = 1
        threshold_count = 0

        qset = StrainMutation.objects.filter(strain_id__in=strains.values('pk')).order_by('mutation_id')

        for pc, mutation_id in chunked_iterator(qset.values_list('mutation_id', flat=True)):

            if st_count % 100000 == 0:
                print(f" [ ] Looping strain mutations {pc:0.2g}% / {st_count} / {mut_count}>{threshold_count}")

            st_count += 1
            if this_mutation == mutation_id:
                this_count += 1
                continue

            if this_mutation is not None:
                bulk_count += 1
                if this_count < threshold:
                    threshold_count += 1
                    continue
                mut_count += 1
                    
                bulk.append(
                    StrainMutationCount(
                        cache_index=index,
                        mutation_id=this_mutation,
                        count=(this_count + 1)
                    )
                )
                this_count = 0
            this_mutation = mutation_id

            if bulk_count > 100:
                sys.stdout.write(f" [x] Bulk creating {pc:0.2g}% ({mut_count} created {threshold_count} mutations ignored)")
                sys.stdout.flush()
                StrainMutationCount.objects.bulk_create(bulk)
                bulk_count = 0
                bulk = []
                sys.stdout.write(f" ... DONE\n")


        if bulk_count:
            print(f" [x] Finishing {mut_count} Mutations")
            StrainMutationCount.objects.bulk_create(bulk)

        return mut_count

    def generate_index(self, index):
        print("Generating index...")
        count, details = index.mutations.all().delete()
        print(f" [x] Deleted {count} existing caching rows.")

        #strains = StrainMutation.objects.filter(index.strain_query('strain__')).order_by()
        #count = strains.count()
        #print(f" [c] Counting {count} StrainMutations into Mutations")

        #strains = strains.filter(mutation=OuterRef('pk')).values('mutation')

        bulk = []
        total_count = 0 
        count = 0 
        #qset_1 = Mutation.objects.annotate(
        #            strain_count=Subquery(strains.annotate(c=Count('*')).values('c'))
        #        ).filter(strain_count__gt=1).values_list('pk', 'strain_count')

        qset = Mutation.objects.filter(index.strain_query('strain_mutations__strain__'))
        with timer() as t:
            qset_2 = qset.annotate(strain_count=Count('strain_mutations__pk'),).filter(strain_count__gt=1)
        print(f"Generate query: {t():.4f} secs")

        #with timer() as t:
        #    print(qset_2.count())
        #print(f"Count rows: {t():.4f} secs")

        with timer() as t:
            print(qset_2.order_by('strain_count')[:100])
        print(f"Get a page of results: {t():.4f} secs")

        sys.exit(1)

        start = perf_counter()
        for x, mutation in enumerate(qset_2[:100]): #chunked_iterator(qset_2, 500):
            if count == 0:
                t = start = perf_counter()
                print(f" [!] First page in {t:.4f}s")
            bulk.append(
                StrainMutationCount(
                    cache_index=index,
                    mutation_id=mutation.pk,
                    count=mutation.strain_count)
                )
            count += 1
            if count >= 100:
                t = start = perf_counter()
                print(f" [x] Created 100 Counts {page:0.2g} in {t:.4f}s")
                StrainMutationCount.objects.bulk_create(bulk)
                total_count += count
                count = 0 
                bulk = []
        if count:
            print(f" [x] Finishing {count} Counts")
            StrainMutationCount.objects.bulk_create(bulk)
        return total_count + count

