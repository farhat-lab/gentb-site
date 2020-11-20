
from django.db import transaction
from django.core.management.commands.loaddata import Command as LoadCommand, parse_apps_and_model_labels
from ...serializers import json_deserializer

class Command(LoadCommand):
    def handle(self, *fixture_labels, **options):
        json_deserializer()

        self.ignore = options['ignore']
        self.using = options['database']
        self.app_label = options['app_label']
        self.verbosity = options['verbosity']
        self.excluded_models, self.excluded_apps = parse_apps_and_model_labels(options['exclude'])
        self.format = options['format']

        self.loaddata(fixture_labels)

        if transaction.get_autocommit(self.using):
            connections[self.using].close()
