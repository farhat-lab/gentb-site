#!/usr/bin/env python3
# pylint: disable=wrong-import-position
# dataset from WHO: LTBI_estimates_2020-07-23.csv

import sys

sys.path.insert(0, '.')
sys.path.insert(0, '..')

try:
    import manage # pylint: disable=unused-import
except ImportError as err:
    sys.stderr.write("Could not run script! Is manage.py not in the current"\
        "working directory, or is the environment not configured?:\n"\
        "{}\n".format(err))
    sys.exit(1)

from apps.mutations.utils import csv_merge
from apps.maps.models import Country, CountryHealth

if __name__ == '__main__':
    with open(sys.argv[1]) as fhl:
        rows = csv_merge(fhl)
        next(rows)
        for row in rows:
            try:
                country = Country.objects.get(iso3=row['iso3'])
                health = country.health
            except CountryHealth.DoesNotExist:
                health = CountryHealth(country=country)
            except Country.DoesNotExist:
                print("Can't find country {}".format(row['iso3']))
                continue
            health.household = row['e_hh_size']
            health.save()
            print("Saved WHO TB data for {}".format(country))
