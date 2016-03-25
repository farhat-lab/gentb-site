"""
This fires off the perl script to kick off the pipeline for either
FastQ or VCF file analysis.

"""
from __future__ import print_function

import logging
from os.path import isfile, join

if __name__ == '__main__':
    import sys
    from os.path import dirname, realpath
    DJANGO_PROJ_DIR = dirname(dirname(dirname(realpath(__file__))))
    sys.path.append(DJANGO_PROJ_DIR)
    #os.environ['DJANGO_SETTINGS_MODULE'] = 'tb_website.settings.local'

    # Allows the working environ to get set-up, apps registered, etc
    #
    import django
    django.setup()

from django.conf import settings
from django.template.loader import render_to_string

from apps.utils.result_file_info import RESULT_FILE_NAME_DICT,\
            EXPECTED_FILE_DESCRIPTIONS,\
            RESULT_OUTPUT_DIRECTORY_NAME
from apps.predict.models import PredictDataset,\
            PredictDatasetNote, DatasetScriptRun,\
            DATASET_STATUS_FILE_RETRIEVAL_COMPLETE

from apps.script_helper.script_runner_basic import run_script

# We could use __file__ here instead
SCRIPT_DIR = join(settings.SITE_ROOT, 'apps', 'predict', 'predict_pipeline')

LOGGER = logging.getLogger('apps.predict.pipeline_hardcoded_script_runner')

class ErrMsg(object):
    """
    Class to hold an error message title and note
    """
    def __init__(self, title, note):
        self.title = title
        self.note = note

VCF_ANALYSIS_SCRIPT = 'analyseVCF.pl'
FASTQ_ANALYSIS_SCRIPT = 'analyseNGS.pl'


class PipelineScriptRunner(object):
    """
    Given a PredictDataset, run the appropriate
    cluster pipeline script to process it.
    """
    def __init__(self, selected_dataset):
        assert isinstance(selected_dataset, PredictDataset),\
            'The dataset must be an instance of PredictDataset'

        self.dataset = selected_dataset

        self.err_found = False
        self.err_message_title = None
        self.err_message = None


    def record_error(self, msg_title, msg):
        """
        Record an error:
            - within this non-persistnt object
            - in the log files
            - in the database in a PredictDatasetNote object
        """
        self.err_found = True
        self.err_message_title = msg_title
        self.err_message = msg

        # Write to the logs
        #
        LOGGER.error(self.err_message_title)
        LOGGER.error(self.err_message)

        # Write to the database
        #
        note = PredictDatasetNote(dataset=self.dataset,\
                      title=self.err_message_title,\
                      note=self.err_message)
        note.save()

        return self.get_err_msg_obj()

    def get_err_msg_obj(self):
        """
        Create an object for sending a message
        back to the calling function
        """
        return ErrMsg(self.err_message_title,\
                    self.err_message)

    def run_script_on_dataset(self):
        """
        Overarching method to excecute:
        - Step 1:  Get script directory information
        - Step 2:  Create/Format the pipeline script command
        - Step 3:  Run the pipeline script!
        """
        script_command = self.get_script_command()
        if script_command is None:
            return (False, self.get_err_msg_obj())

        return self.run_command(script_command)


    def get_script_command(self):
        """
        Using dataset information, to decide whether to run:
            (1) script for a VCF file
            (2) script for FastQ files
        """
        # Formate either a VCF or FastQ pipeline command
        #
        if self.dataset.is_vcf_file():     # (2a) command for a VCF file
            command_to_run = self.get_vcf_script_command()

        elif self.dataset.is_fastq_file(): # (2b) command for FastQ files
            command_to_run = self.get_fastq_script_command()

        else:
            err_title = 'Not VCF or FastQ file'
            err_note = 'Could not determine the file type.\
             Database contained: "%s"' % (self.dataset.file_type)
            #err_msg_obj =
            self.record_error(err_title, err_note)
            return None

        return command_to_run



    def run_command(self, command_to_run):
        """
        - (1) Create a DatasetScriptRun object
        - (2) Make a custom callback script and put it in the dataset file directory
        - (3) Update the Dataset status
        - (4) Run the script
        """
        if command_to_run is None:
            return (False, self.get_err_msg_obj())

        # (1) Create a run object and save the command being run
        #
        dsr = DatasetScriptRun(dataset=self.dataset)
        dsr.notes = command_to_run
        dsr.save()

        # ---------------------------------------------------------
        # (2) Place a callback script with the Dataset directory
        #   - contains callback url + md5
        #   - called after the analyse scripts are Run
        # ---------------------------------------------------------

        # Get callback args
        callback_info_dict = self.dataset.get_script_args_json(dsr.md5, as_dict=True)
        template_dict = dict(callback_info_dict=callback_info_dict)
        template_dict.update(RESULT_FILE_NAME_DICT)
        template_dict['RESULT_OUTPUT_DIRECTORY_NAME'] = RESULT_OUTPUT_DIRECTORY_NAME
        template_dict['EXPECTED_FILE_DESCRIPTIONS'] = EXPECTED_FILE_DESCRIPTIONS

        print('-' * 40)
        print(template_dict)
        print('-' * 40)

        # Use args to create a custom python script
        py_callback_script = render_to_string('feedback/gentb_status_feedback.py', template_dict)

        # Place script in dataset file_directory
        script_fullname = join(self.dataset.file_directory, 'gentb_status_feedback.py')
        fhandler = open(script_fullname, 'w')
        fhandler.write(py_callback_script)
        fhandler.close()

        # (3) Update the dataset status to 'in process'
        #
        self.dataset.set_status_processing_started()

        # (4) Run the script -- not waiting for output
        #
        full_args = command_to_run.split()
        print ('full_args', full_args)
        run_script(full_args)

        return (True, dsr)


    def get_fastq_script_command(self):
        """
        Run the command for FastQ files.
        e.g. perl analyseNGS.pl (directory name)
        """
        if self.err_found:
            return None

        # (1) Make sure the 'analyseNGS.pl' command is in the
        #   specified 'Pipeline Scripts Directory'
        #
        script_cmd = join(SCRIPT_DIR, FASTQ_ANALYSIS_SCRIPT)
        if not isfile(script_cmd):
            err_title = 'The "%s" file was not found' % (FASTQ_ANALYSIS_SCRIPT)
            err_note = """The "%s" file was not found at this location: \n%s\n
            To change the path to the script, please go into the admin control\
             panel and modify the 'Pipeline Scripts Directory'.\n\
            Talk to your administrator for details.""" % (FASTQ_ANALYSIS_SCRIPT, script_cmd)
            #err_msg_obj = self.record_error(err_title, err_note)
            self.record_error(err_title, err_note)
            return None

        # Format the full command with target containing
        #   input files
        #
        if self.dataset.is_fastq_single_ended():
            command_str = 'perl {0} 0 . {1}'.format(script_cmd,\
                self.dataset.file_directory)
        elif self.dataset.is_fastq_pair_ended():
            pair_extension = self.dataset.get_fastq_pair_end_extension()
            if pair_extension is None:
                err_title = 'FastQ could not find pair-ended extension type'
                err_note = 'Could not determine pair-ended extension type.\
                 Database contained: "%s"' % (self.dataset.fastq_type)
                #err_msg_obj = self.record_error(err_title, err_note)
                self.record_error(err_title, err_note)
                return None

            command_str = 'perl {0} 1 {1} {2}'.format(script_cmd,\
                pair_extension,\
                self.dataset.file_directory)
        else:
            err_title = 'FastQ: single-ended or pair-ended?'
            err_note = 'Could not determine single-ended or pair-ended FastQ type.\
             Database contained: "%s"' % (self.dataset.fastq_type)
            # err_msg_obj = self.record_error(err_title, err_note)
            self.record_error(err_title, err_note)
            return None

        return command_str



    def get_vcf_script_command(self):
        """
        Run the command for a VCF file.
        e.g. perl analyseVCF.pl (directory name)
        """
        if self.err_found:
            return None

        # (1) Make sure the 'analyseVCF.pl' command is in the
        #   specified 'Pipeline Scripts Directory'
        #
        script_cmd = join(SCRIPT_DIR, VCF_ANALYSIS_SCRIPT)
        if not isfile(script_cmd):
            err_title = 'The "%s" file was not found' % (VCF_ANALYSIS_SCRIPT)
            err_note = """The "%s" file was not found at this location: \n%s\n
            To change the path to the script, please go into the admin control\
             panel and modify the 'Pipeline Scripts Directory'.\n\
            Talk to your administrator for details.""" % (VCF_ANALYSIS_SCRIPT, script_cmd)
            #err_msg_obj = self.record_error(err_title, err_note)
            self.record_error(err_title, err_note)
            return None

        # Format the full command with target containing
        #   input files
        #
        command_str = 'perl {0} {1}'.format(script_cmd,\
                self.dataset.file_directory)

        return command_str

    @staticmethod
    def get_pipeline_command(selected_dataset):
        """
        Return the full pipeline command to run the
        appropriate analyze script for this dataset's files
        """
        if selected_dataset is None:
            return (False, 'The selected_dataset is None')

        pipeline_runner = PipelineScriptRunner(selected_dataset)

        script_command = pipeline_runner.get_script_command()
        if script_command is None:
            return (False, pipeline_runner.err_message)

        return (True, script_command)

    @staticmethod
    def run_next_dataset():
        """
        (1) Retrieve the first PredictDataset with a status of:
            DATASET_STATUS_FILE_RETRIEVAL_COMPLETE
        (2) If such a dataset exists, run it through the pipeline
        """
        #pass
        LOGGER.debug("Run pipeline check: next dataset")
        print("Run pipeline check: next dataset")

        # get some Dataset
        dataset_from_db = PredictDataset.objects.filter(\
                status=DATASET_STATUS_FILE_RETRIEVAL_COMPLETE).first()

        if dataset_from_db is None:
            LOGGER.debug("Nothing to check")
            print('Nothing to check')
            return

        # Run script
        LOGGER.debug("Run pipeline for dataset: %s (%s)",\
            dataset_from_db, dataset_from_db.id)
        print("Run pipeline for dataset: %s (%s)" %\
            (dataset_from_db, dataset_from_db.id))
        pipeline_runner = PipelineScriptRunner(dataset_from_db)
        pipeline_runner.run_script_on_dataset()


if __name__ == '__main__':
    PipelineScriptRunner.run_next_dataset()

    """
    # Alternative run method
    script_command = pipeline_runner.get_script_command()
    if script_command is not NOne:
        pipeline_runner.run_command(script_command)


    """
