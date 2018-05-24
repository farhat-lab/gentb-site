#
# Copyright (C) 2017 Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Test the dropbox upload functionality.

We can't really test the javascript from here, so we'll focus
on making sure the backend get-files portion works.
"""
import os
import json

from django.contrib.auth import get_user_model

from autotest.base import ExtraTestCase

from ..forms import TestUploadForm
from ..models import DropboxUploadFile

class TestDropbox(ExtraTestCase):

    def test_upload(self):
        fn = os.path.join(self.app_dir, 'tests', 'media', 'todo.vcf')
        self.assertUpload('user217', 'dropbox', link='file://' + fn)

    def assertUpload(self, username, source, **attr):
        get_user_model().objects.create(username='r217')
        self.client.login(username='r217', password=True)

        data = {
          "id":"todo.vcf",
          "name":"todo.vcf",
          "bytes":300,
          "source": source,
          "bucket": None,
          "icon": "test_icon.png",
          }
        data.update(attr)

        response = self.assertPost('uploads:test', data=dict(test_files=json.dumps([data])))

        # This depends on the form existing in the response
        form = self.contextData(response, 'form')
        self.assertTrue(form.is_valid())

        files = form.cleaned_data.get('test_files')
        self.assertEqual(len(files[None]), 1)
        self.assertEqual(type(files[None][0]), DropboxUploadFile)
        self.assertEqual(files[None][0].pk, None)

        f = files[None][0]
        f.conclude_upload(self.media_root, self.user)
        self.assertTrue(f.pk)

        self.assertEqual(f.retrieval_start, None)
        self.assertEqual(f.retrieval_end, None)
        self.assertEqual(f.retrieval_error, '')

        fn = os.path.join(self.media_root, f.filename)
        self.assertFalse(os.path.isfile(fn))

        f.download_now()

        self.assertEqual(f.retrieval_error, '')
        self.assertTrue(f.retrieval_start)
        self.assertTrue(f.retrieval_end)
        self.assertTrue(os.path.isfile(fn))

    def test_retry_upload(self):
        pass
