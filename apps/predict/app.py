#
# Copyright 2018, Maha Farhat
#
# This file is part of the software gentb, consisting of custom
# code for the Inkscape project's django-based website.
#
# gentb is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gentb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with gentb.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Signals for pipeline
"""

from django.db.models import signals
from django.apps import AppConfig

class PredictApp(AppConfig):
    """Predict Application with signals"""
    name = 'apps.predict'

    @staticmethod
    def delete_piperun(instance, **_):
        """Remove any pipelines associated with this predict strain"""
        # Pipeline runs should be deleted with the predict dataset.
        if instance.piperun:
            instance.piperun.delete()
        for upload in instance.files:
            upload.delete()

    def ready(self):
        """Called when the app is ready"""
        from .models import PredictStrain
        signals.pre_delete.connect(self.delete_piperun, sender=PredictStrain)
