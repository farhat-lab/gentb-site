"""
Load a custom map via csv
"""

import os
import csv
import json
import sys

from django.core.management.base import BaseCommand, CommandError

from apps.maps.utils import COUNTRY_MAP
from apps.maps.models import Country
from apps.mutations.models import Drug
from apps.mutations.utils import csv_merge
from apps.maptables.models import CustomMap, MapRow

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('filename', type=str)
        parser.add_argument('key', type=str)
        parser.add_argument('--keep', action='store_true')

    def get_delim(self, filename):
        if filename.endswith('csv'):
            return ','
        if filename.endswith('tsv'):
            return '\t'
        raise CommandError("Unknown csv type (expected csv or tsv)")

    def handle(self, filename, key, keep=False, **options):
        if not os.path.isfile(filename):
            raise CommandError(f"File not found: {filename}")
        delim = self.get_delim(filename)

        custom_map, created = CustomMap.objects.get_or_create(
            slug=key,
            defaults={'name': key, 'description': 'Automatically created'}
        )
        if created:
            print(f"Created new custom map {key}...")
        elif keep:
            print(f"Adding to existing map {key}...")
        else:
            rows, _ = custom_map.rows.all().delete()
            print(f"Deleted {rows} rows in {key} and reloading...")

        rows_added = 0
        bad_drugs = set()
        bad_countries = set()
        with open(filename, 'r') as fhl:
            for row in csv_merge(fhl, delimiter=delim):
                if isinstance(row, list):
                    if 'drug' not in row:
                        raise CommandError("Can't find column 'drug' in csv.")
                    if 'country' not in row:
                        raise CommandError("Can't find column 'country' in csv.")
                elif isinstance(row, dict):
                    try:
                        self.add_row(custom_map, row)
                        rows_added += 1
                    except Country.DoesNotExist:
                        bad_countries.add(row['country'])
                    except Drug.DoesNotExist:
                        bad_drugs.add(row['drug'])

        if bad_drugs:
            print(f" * WARNING: Could not match drugs {bad_drugs}, rows not imported")
        if bad_countries:
            print(f" * WARNING: Could not match countries {bad_countries}, rows not imported")
        if not rows_added:
            print(f" * CRITICAL: No rows added!")
        else:
            print(f"\n = Sucessfully imported {rows_added} rows =\n")


    def add_row(self, custom_map, row):

        drug = Drug.objects.get(code__iexact=row['drug'].lower())

        # Try everything to get the country
        c_name = COUNTRY_MAP.get(row['country'], row['country'])
        country = None
        for c_col in ('iso2', 'iso3', 'name__iexact'):
            try:
                country = Country.objects.get(**{c_col: c_name})
                break
            except Country.DoesNotExist:
                continue
        if country is None:
            raise Country.DoesNotExist("Can't find country")

        cleaned = {}
        for col, val in row.items():
            if col in ('drug', 'country'):
                continue

            if val.isnumeric():
                cleaned[col] = int(val)
            elif val.count('.') == 1 and val.replace('.', '').isnumeric():
                cleaned[col] = float(val)
            else:
                cleaned[col] = row[col]

        custom_map.rows.update_or_create(
            drug=drug,
            country=country,
            defaults={'data': json.dumps(cleaned)})

