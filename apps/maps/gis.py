"""
Provide a little bit of abstraction for gis.
"""
try:
    from django.contrib.gis.db.models import *
except Exception:
    def MultiPolygonField(*args, **kw):
        return TextField()
    MultiPointField = MultiPolygonField
    GeoManager = Manager
