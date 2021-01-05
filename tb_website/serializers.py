"""
Serialise very large data tables for import.
"""

import os
import sys
import json
import atexit

from collections import defaultdict

from django.db import models, transaction
from django.core.serializers import base, get_serializer_formats, register_serializer
from django.core.serializers.python import _get_model
from django.core.serializers.json import Serializer as BaseSerializer
from django.core.exceptions import ObjectDoesNotExist

from .utils import sizeof, to

class Serializer(BaseSerializer):
    pass

class ProgressiveLoader():
    reject_fhl = None

    def __init__(self, fhl, pos=1):
        self.disk = fhl
        self.pos = pos

    def __iter__(self):
        """Loop through disk data"""
        self.disk.seek(self.pos) # Skip '[' == 1
        buf = b''
        while True:
            self.pos = self.disk.tell()
            chunk = self.disk.read(200000)
            buf += chunk

            while b"\n}" in buf:
                (a, buf) = buf.split(b"\n{\n", 1)
                (segment, buf) = buf.split(b"\n}", 1)
                buf = buf.lstrip(b",")
                self.block = json.loads("{\n" + segment.decode('utf8') + "\n}")
                yield self.block['model'], self.block.get('pk', None), self.block['fields']
                self.pos = self.disk.tell() - len(buf)

            if not chunk:
                break

    def reject_current_block(self):
        if not self.block:
            return
        self.reject(json.dumps(self.block, indent=2), self.disk.name + '.rejected')

    @classmethod
    def reject(cls, block, base_name):
        prefix = ""
        if not cls.reject_fhl:
            count = 0
            filename = base_name
            while os.path.isfile(filename):
                filename = base_name + str(count)
                count += 1
            cls.reject_fhl = open(filename, 'w')
            cls.reject_fhl.write("[\n")
            atexit.register(cls.close_reject)
        else:
            prefix = ",\n"
        cls.reject_fhl.write(prefix + block + "\n")

    @classmethod
    def close_reject(cls):
        if cls.reject_fhl:
            cls.reject_fhl.write("]")
            cls.reject_fhl.close()

class BigDeserializer():
    field_names_cache = {}

    def __init__(self, ranges=None, attempt=1):
        self.has_ranges = bool(ranges)
        self.ranges = ranges
        self.attempt = attempt
        # A list of models which have had problems loaded and are now being deferred,
        # this uses a huristic weighting where errors increment and successes decrement.
        self.deferred_models = defaultdict(int)

        # Where we have items which are deferred, we want to record their indexes.
        self.deferred_ranges = []

    def __call__(self, fhl, pos=0, **options):
        """We are replacing the json deserializer completly."""
        if isinstance(fhl, (bytes, str)):
            raise IOError("String or bytes not accepted, only file handles")

        count = 0
        if self.has_ranges and pos != 0:
            count = self.ranges[0][0] - 1
        success = 0
        last_pos = 0
        problem_pos = -1
        size = os.path.getsize(fhl.name)
        loader = ProgressiveLoader(fhl, pos=pos)
        for m_name, pk, fields in loader:
            count += 1
            if fhl.tell() != last_pos:
                sys.stdout.write("Loading {}: {}/{} ({:.2f}%) {} done, {} deferred\r".format(
                    self.attempt, sizeof(loader.pos), sizeof(size),
                    loader.pos / size * 100, success, count - success))
                last_pos = fhl.tell()

            while self.ranges and count > self.ranges[0][1]:
                self.ranges.pop(0)

            if (self.has_ranges and not self.ranges) or (self.ranges and count < self.ranges[0][0]):
                success += 1
                continue # Out of range item

            obj = self.load_block(m_name, pk, fields)
            if obj is not None:
                try:
                    obj.save()
                except Exception:
                    # Anything that can go wrong.
                    obj = None

            if obj is not None:
                #yield obj
                success += 1
                if success % 1000 == 0:
                     transaction.commit()
            else:
                if not self.has_ranges:
                    loader.reject_current_block()
                problem_pos = loader.pos
                if self.deferred_ranges and self.deferred_ranges[-1][1] == count - 1:
                    self.deferred_ranges[-1][1] = count
                else:
                    self.deferred_ranges.append([count, count])

        print("\n")
        # Some successes but also problems, retry.
        if success > 0 and problem_pos != -1:
            yield from BigDeserializer(self.deferred_ranges, self.attempt+1)(fhl, pos=problem_pos)



    def load_block(self, m_name, pk, data):
        """
        Import a single block of fields into the given model.
        """
        if self.deferred_models[m_name] > 3:
            # Model is deferred!
            return None

        try:
            model, fields = self.get_model(m_name)
        except base.DeserializationError:
            sys.stderr.write("\nSkipping all items from model {m_name}!\n")
            self.deferred_models[m_name] = 10
            return None

        try:
            data = self.build_data(model, fields, pk, data)
            m2m_fields = data.pop('m2m', {})
        except (base.M2MDeserializationError, base.DeserializationError, ObjectDoesNotExist):
            return None

        try:
            obj = base.build_instance(model, data, None)
            return base.DeserializedObject(obj, m2m_fields, {})
        except Exception:
            return None # Any error, defer the object.

    @to(dict)
    @staticmethod
    def build_data(model, field_names, pk, data, using=None):
        m2m = {}
        if pk is not None:
            try:
                yield (model._meta.pk.attname, model._meta.pk.to_python(pk))
            except Exception as err:
                raise base.DeserializationError.WithData(err, m_name, pk, None)

        for (field_name, field_value) in data.items():
            if field_name not in field_names:
                continue

            field = model._meta.get_field(field_name)

            if field.remote_field and isinstance(field.remote_field, models.ManyToManyRel):
                values = base.deserialize_m2m_values(field, field_value, using, False)
                m2m[field.name] = values

            # Handle FK fields
            elif field.remote_field and isinstance(field.remote_field, models.ManyToOneRel):
                value = base.deserialize_fk_value(field, field_value, using, False)
                yield (field.attname, value)

            else:
                # Everything else
                yield (field.name, field.to_python(field_value))
        # Not great
        yield ('m2m', m2m)

    def get_model(self, name):
        """Gets the model name and fields"""
        model = _get_model(name)
        if model not in self.field_names_cache:
            self.field_names_cache[model] = {f.name for f in model._meta.get_fields()}
        return (model, self.field_names_cache[model])

def json_deserializer():
    get_serializer_formats()
    from django.core.serializers import _serializers
    register_serializer('json', 'tb_website.serializers', _serializers)

Deserializer = BigDeserializer()
