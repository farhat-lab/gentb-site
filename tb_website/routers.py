
from django.conf import settings
DB = list(settings.DATABASES)

class PrivateAppRouter():
    """
    Allows an app to have it's own private database if needed
    """
    @staticmethod
    def db_for_read(model, **hints):
        if model._meta.app_label in DB:
            return model._meta.app_label
    db_for_write = db_for_read

    @staticmethod
    def allow_relation(obj1, obj2, **hints):
        if obj1._state.db == obj2._state.db:
            return True
        return None

    @staticmethod
    def allow_migrate(db, app_label, model_name=None, **hints):
        return app_label not in DB

