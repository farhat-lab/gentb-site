from os.path import dirname, realpath
import json

from django.test import TestCase
from apps.dropbox_helper.forms import DropboxRetrievalParamsForm


class RetrievalParamsTestCase(TestCase):
    def setUp(self):
        self.test_params = dict(\
            dropbox_url='https://a-dropbox-shared-link.com',
            destination_directory=dirname(realpath(__file__)), callback_url='https://myserver.com/predict/file-retrieval-results', callback_md5='8fec9fefa93095fc94a68f495e24325b'\
            )


    def test_params_form(self):

        #print json.dumps(self.test_params)
        # -----------------
        # Good params
        # -----------------
        f = DropboxRetrievalParamsForm(self.test_params)
        self.assertTrue(f.is_valid(),
                        "Form failed: {0}".format(f.errors))

        # -----------------
        # Bad params 1
        # -----------------
        bad_params = self.test_params.copy()
        bad_params['callback_md5'] = 'too-short'

        f = DropboxRetrievalParamsForm(bad_params)
        self.assertTrue(f.is_valid() is False)
        self.assertTrue('callback_md5' in f.errors.keys())
        md5_err = [u'Not a valid MD5. Must be 32 characters.']
        self.assertEqual(f.errors.values()[0], md5_err )

        # -----------------
        # Bad params 2
        # -----------------
        bad_params2 = self.test_params.copy()
        bad_params2['destination_directory'] = '/no-dir-here/we-hope'

        f = DropboxRetrievalParamsForm(bad_params2)
        self.assertTrue(f.is_valid() is False)


        self.assertTrue('destination_directory' in f.errors.keys())
        dir_err = [u'Destination directory not found: /no-dir-here/we-hope']
        self.assertEqual(f.errors.values()[0], dir_err )
