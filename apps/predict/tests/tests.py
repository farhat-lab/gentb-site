from django.test import TestCase

from os.path import dirname, realpath, join
import json

from django.contrib.auth.models import User
from django.utils.crypto import get_random_string

from apps.predict.models import PredictDataset
from apps.tb_users.models import TBUser
from apps.utils.file_patterns import \
                                GENTB_FASTQ_FILE_PATTERNS,\
                                GENTB_VCF_FILE_PATTERNS,\
                                FILE_TYPE_VCF,\
                                FILE_TYPE_FASTQ,\
                                FASTQ_PAIR_ENDED


class PredictBasicTest(TestCase):
    """
    Run basic tests with the PredictDataset object
    """
    def setUp(self):
        """
        Create initial objects for testing
        """
        super(PredictBasicTest, self).setUp()
        self.test_user = User(username='testtb', password=get_random_string(length=32))
        self.test_user.save()

        test_params = {
            u'description': u'ok',\
            u'file_type': u'vcf',\
            u'title': u'vcf - embed',\
            u'file_directory':\
                u'/some-dir-to-add-files/test_setup/tb_uploaded_files/tbdata_00000013',\
            u'user': self.test_user}

        self.dataset_vcf = PredictDataset(**test_params)
        self.dataset_vcf.save()

        test_params2 = test_params.copy()
        test_params2['file_type'] = FILE_TYPE_FASTQ
        self.dataset_fastq = PredictDataset(**test_params2)
        self.dataset_fastq.save()

    def test_prediction_pipeline(self):
        pass

