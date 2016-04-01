from django.test import TestCase

from os.path import dirname, realpath, join
import json

from django.contrib.auth.models import User
from django.utils.crypto import get_random_string

from apps.predict.models import PredictDataset, PredictDatasetStatus,\
    SCRIPT_DIR, DATASET_STATUS_FILE_RETRIEVAL_COMPLETE
from apps.tb_users.models import TBUser
from apps.utils.file_patterns import FilePatternHelper,\
                                GENTB_FASTQ_FILE_PATTERNS,\
                                GENTB_VCF_FILE_PATTERNS,\
                                FILE_TYPE_VCF,\
                                FILE_TYPE_FASTQ,\
                                FASTQ_PAIR_ENDED


class PredictBasicTest(TestCase):
    """
    Run basic tests with the PredictDataset object
    """
    fixtures = ['initial_data']

    def setUp(self):
        """
        Create initial objects for testing
        """
        self.test_user = User(username='testtb', password=get_random_string(length=32))
        self.test_user.save()

        self.tb_test_user = TBUser(user=self.test_user, affiliation='HU')
        self.tb_test_user.save()

        test_params = {u'status':\
                PredictDatasetStatus.objects.get(pk=DATASET_STATUS_FILE_RETRIEVAL_COMPLETE),\
            u'has_prediction': False,\
            u'description': u'ok',\
            u'file_type': u'vcf',\
            u'title': u'vcf - embed',\
            u'file_directory':\
                u'/some-dir-to-add-files/test_setup/tb_uploaded_files/tbdata_00000013',\
            u'fastq_type': u'',\
            u'dropbox_url': \
                u'https://www.dropbox.com/sh/p6ses8376312bes/AAA7TB4GhErfLLfE7WPco79ha?dl=0',\
            u'user': self.tb_test_user}

        self.dataset_vcf = PredictDataset(**test_params)
        self.dataset_vcf.save()

        test_params2 = test_params.copy()
        test_params2['file_type'] = FILE_TYPE_FASTQ
        test_params2['fastq_type'] = FASTQ_PAIR_ENDED
        self.dataset_fastq = PredictDataset(**test_params2)
        self.dataset_fastq.save()

    def tearDown(self):
        """
        Delete test objects
        """
        self.dataset_vcf.delete()
        self.dataset_fastq.delete()
        self.tb_test_user.delete()
        self.test_user.delete()
        self.pipleline_scripts_info.delete()

    def test_vcf_file(self):
        """VCF file test"""
        self.assertTrue(self.dataset_vcf.is_vcf_file(),
                        "Should be a VCF file")

        self.assertTrue(not self.dataset_vcf.is_fastq_file(),
                        "Should not be FastQ file")

        self.assertEqual(self.dataset_vcf.get_file_patterns(),\
                    GENTB_VCF_FILE_PATTERNS)

        self.assertTrue(self.dataset_vcf.get_pipeline_command()[0])

        self.dataset_vcf.run_command()

    def test_fastq_file(self):
        """FastQ file test"""
        self.assertTrue(not self.dataset_fastq.is_vcf_file(),
                        "Should NOT be a VCF file")

        self.assertTrue(self.dataset_fastq.is_fastq_file(),
                        "Should be a FastQ file")

        self.assertEqual(self.dataset_fastq.get_file_patterns(),\
                    GENTB_FASTQ_FILE_PATTERNS)

        self.assertTrue(self.dataset_fastq.get_pipeline_command()[0])

        self.dataset_fastq.run_command()

