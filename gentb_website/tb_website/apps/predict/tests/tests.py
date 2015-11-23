from django.test import TestCase

from os.path import dirname, realpath, join
import json

from django.contrib.auth.models import User
from django.utils.crypto import get_random_string

from apps.predict.models import PredictDataset, PredictDatasetStatus,\
    PipelineScriptsDirectory,\
    DATASET_STATUS_FILE_RETRIEVAL_COMPLETE
from apps.tb_users.models import TBUser
from apps.predict.pipeline_hardcoded_script_runner import PipelineScriptRunner
from apps.utils.file_patterns import FilePatternHelper,\
                                GENTB_FASTQ_FILE_PATTERNS,\
                                GENTB_VCF_FILE_PATTERNS


class PredictBasicTest(TestCase):

    fixtures = ['initial_data']

    def setUp(self):

        self.test_user = User(username='testtb', password=get_random_string(length=32))
        self.test_user.save()

        self.tb_test_user = TBUser(user=self.test_user, affiliation='HU')
        self.tb_test_user.save()

        test_params = {u'status': PredictDatasetStatus.objects.get(pk=DATASET_STATUS_FILE_RETRIEVAL_COMPLETE),
            u'has_prediction': False,
            u'description': u'ok',
            u'file_type': u'vcf',
            u'title': u'vcf - embed',
            u'file_directory': u'/some-dir-to-add-files/test_setup/tb_uploaded_files/tbdata_00000013', u'fastq_type': u'',
            u'dropbox_url': u'https://www.dropbox.com/sh/p6ses8376312bes/AAA7TB4GhErfLLfE7WPco79ha?dl=0',
            u'user': self.tb_test_user }

        self.dataset1 = PredictDataset(**test_params)
        self.dataset1.save()

        script_dirname = join(dirname(realpath(__file__)), 'test_script_dir')
        self.pipleline_scripts_info = PipelineScriptsDirectory(name='test dir',\
                    script_directory=script_dirname)
        self.pipleline_scripts_info.save()

    def tearDown(self):
        self.dataset1.delete()
        self.tb_test_user.delete()
        self.test_user.delete()
        self.pipleline_scripts_info.delete()

    def test_params_form(self):


        self.assertTrue(self.dataset1.is_vcf_file(),
                        "Should be a VCF file")

        self.assertTrue(not self.dataset1.is_fastq_file(),
                        "Should not be FastQ file")

        self.assertEqual(self.dataset1.get_file_patterns(),\
                    GENTB_VCF_FILE_PATTERNS)


        pipeline_runner = PipelineScriptRunner(self.dataset1)

        script_directory = pipeline_runner.step1_get_script_directory_info()
        self.assertTrue(script_directory is not None)
        script_command = pipeline_runner.step2_get_script_command(script_directory)
        print script_command
        self.assertTrue(script_command is not None)
        #if script_command is not NOne:
        #    pipeline_runner.step3_run_command(script_command)



        # Add tests to get the bsub command
