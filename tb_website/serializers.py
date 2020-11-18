"""
Serialise very large data tables for import.
"""

import os
import sys
import json

from django.core.serializers.base import DeserializationError
from django.core.serializers import get_serializer_formats, register_serializer
from django.core.serializers.json import (
  Serializer as JsonSerializer,
  PythonDeserializer,
)

class Serializer(JsonSerializer):
    """No changes to output"""
    pass

class ProgressiveLoader(object):
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
                yield json.loads("{\n" + segment.decode('utf8') + "\n}")
                self.pos = self.disk.tell() - len(buf)

            if not chunk:
                break

from math import log
unit_list = ['bytes', 'kB', 'MB', 'GB', 'TB', 'PB']
decs = [0, 0, 1, 2, 2, 2]
def sizeof_fmt(num):
    """Human friendly file size"""
    if num > 1:
        exponent = min(int(log(num, 1024)), len(unit_list) - 1)
        quotient = float(num) / 1024**exponent
        unit = unit_list[exponent]
        num_decimals = decs[exponent]
        format_string = '{:.%sf} {}' % (num_decimals)
        return format_string.format(quotient, unit)
    if num == 0:
        return '0 bytes'
    if num == 1:
        return '1 byte'

def Deserializer(fhl, **options):
    """We are replacing the json deserializer completly."""
    if isinstance(fhl, (bytes, str)):
        raise IOError("String or bytes not accepted, only file handles")

    try:
        count = 0
        last_pos = 0
        size = os.path.getsize(fhl.name)
        loader = ProgressiveLoader(fhl)
        descer = PythonDeserializer(loader, **options)
        while True:
            count += 1
            try:
                obj = next(descer)
            except StopIteration:
                print(f"Complete {count}")
                break
            except Exception as err:
                print(f"\no:{count} ! Exception: {err}")
                descer = PythonDeserializer(loader, **options)
                continue
            if fhl.tell() != last_pos:
                sys.stdout.write("Loading: {}/{} ({:.2f}%) {} objects\r".format(
                    sizeof_fmt(loader.pos), sizeof_fmt(size), loader.pos / size * 100, count))
                last_pos = fhl.tell()
            yield obj
    except GeneratorExit:
        print("\n\n")
        raise
    except Exception as e:
        print("\n\n")
        raise DeserializationError(str(e))

def json_deserializer():
    get_serializer_formats()
    from django.core.serializers import _serializers
    register_serializer('json', 'tb_website.serializers', _serializers)

