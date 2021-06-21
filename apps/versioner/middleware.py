#
# Copyright 2015, Martin Owens
#           2017, MIT TESS
#           2017, Maha Farhat
#
# versioner is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# versioner is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with versioner.  If not, see <http://www.gnu.org/licenses/>.
#
# Code imported from inkscape-web 2017-10-18 AGPLv3
#

try:
    from git import Repo
except (ImportError, ModuleNotFoundError):
    Repo = None

from django.conf import settings
from .utils import BaseMiddleware, to

BRANCHES = getattr(settings, 'VERSION_BRANCHES', [])
ORIGIN = getattr(settings, 'VERSION_ORIGIN', 'origin')

# XXX Future feature, auto-fetch the data when needed, so we
# don't need a cron job.
#FETCH = getattr(settings, 'VERSION_FETCH_EVERY', None)

class VersionInformation(BaseMiddleware):
    """ 
    Load any version information cached on the disk.
    """
    def process_template_response(self, request, response):
        if not hasattr(response, 'context_data') or Repo is None:
            return response

        self.repo = Repo(settings.SITE_ROOT)
        local = self.repo.active_branch.name

        branches = []
        for label, remote in BRANCHES:
            branch = self.get_iter(local, remote)
            branches.append((label, branch))
            if branch:
                response.context_data['branch_updates'] = True
    
        response.context_data['branches'] = branches
        return response

    @to(list)
    def get_iter(self, local, remote):
        """
        Gets a list of commits which have changed between local and remote.
        """
        if '/' not in remote:
            remote = '/'.join([ORIGIN, remote])
        git_range = '{local}..{remote}'.format(local=local, remote=remote)
        return self.repo.iter_commits(git_range)

