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

class PipelineApp(AppConfig):
    """Pipeline Application with signals"""
    name = 'apps.pipeline'

    @staticmethod
    def clean_run(instance, **_):
        """Remove any files in an program run instance when deleting"""
        instance.delete_output_files()

    @staticmethod
    def clean_pipe(instance, **_):
        """Remove any files in an program run instance when deleting (make sure!)"""
        for run in instance.programs.all():
            run.delete_output_files()

    def ready(self):
        """Called when the app is ready"""
        from .models import ProgramRun, PipelineRun
        signals.pre_delete.connect(self.clean_pipe, sender=PipelineRun)
        signals.pre_delete.connect(self.clean_run, sender=ProgramRun)
