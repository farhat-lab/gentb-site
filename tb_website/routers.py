
from django.conf import settings
DB = list(settings.DATABASES)

class PrivateAppRouter(object):
    """
    Allows an app to have it's own private database if needed
    """
    def db_for_read(self, model, **hints):
        if model._meta.app_label in DB:
            return model._meta.app_label
    db_for_write = db_for_read

    def allow_relation(self, obj1, obj2, **hints):
        if obj1._state.db == obj2._state.db:
            return True
        return None

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        return app_label not in DB

