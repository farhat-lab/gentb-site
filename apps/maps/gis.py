"""
Provide a little bit of abstraction for gis.
"""
from django.db.models import TextField, Manager

try:
    from django.contrib.gis.db.models import MultiPolygonField, MultiPointField
except Exception:
    raise
    def MultiPolygonField(*args, **kw):
        return TextField()
    MultiPointField = MultiPolygonField

GeoManager = Manager
