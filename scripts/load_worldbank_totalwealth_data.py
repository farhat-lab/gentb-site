#!/usr/bin/env python3
# pylint: disable=wrong-import-position
# dataset using total wealth variable from world bank wealth accounting dataset

import sys
import csv

sys.path.insert(0, '.')
sys.path.insert(0, '..')

try:
    import manage # pylint: disable=unused-import
except ImportError as err:
    sys.stderr.write("Could not run script! Is manage.py not in the current"\
        "working directory, or is the environment not configured?:\n"\
        "{}\n".format(err))
    sys.exit(1)
from apps.maps.models import Country, CountryHealth

if __name__ == '__main__':
    with open(sys.argv[1]) as fhl:
        rows = csv.reader(fhl)
        header = next(rows)
        country_index = header.index("Country Code")
        year_index = header.index("2014 [YR2014]")
        for row in rows:
            try:
                country = Country.objects.get(iso3=row[country_index])
                health = country.health
            except CountryHealth.DoesNotExist:
                health = CountryHealth(country=country)
            except Country.DoesNotExist:
                print("Can't find country {}".format(row[country_index]))
                continue
            health.total_wealth = row[year_index]
            health.save()
            print("Saved World Bank data for {}".format(country))
