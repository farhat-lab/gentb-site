
from django.core.management.commands.loaddata import Command as LoadCommand
from ...serializers import json_deserializer

class Command(LoadCommand):
    def handle(self, *args, **options):
        json_deserializer()
        super(Command, self).handle(*args, **options)

