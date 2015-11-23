"""
This fires off the "bsub" command to kick off the pipeline for either
FastQ or VCF file analysis.

"""
from apps.predict.models import PredictDatasetStatus, PredictDataset,\
            PredictDatasetNote, DatasetScriptRun,\
            PipelineScriptsDirectory
from apps.script_helper.script_runner_basic import run_script
from os.path import isfile, isdir, join
import logging
LOGGER = logging.getLogger('apps.predict.pipeline_hardcoded_script_runner')

class ErrMsg:
    def __init__(self, title, note):
        self.title = title
        self.note = note

class PipelineScriptRunner(object):

    def __init__(self, dataset):
        assert isintance(dataset, PredictDataset),\
            'The dataset must be an instance of PredictDataset'

        self.dataset = dataset

        self.err_found = False
        self.err_message_title = None
        self.err_message = None


    def record_error(self, msg_title, msg):
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
                      title=err_title,\
                      note=err_note)
        note.save()

        return self.get_err_msg_obj()

    def get_err_msg_obj(self):
        # Create an object for sending a message
        # back to the calling function
        #
        return ErrMsg(self.err_message_title,\
                    self.err_message)

    def get_script_directory_info(self):
        """
        Retrieve PipelineScriptsDirectory object from the database
        """
        script_directory_info = PipelineScriptsDirectory.objects.filter(is_chosen_script=True).first()
        if script_directory_info is None:
            err_title='No pipeline directory specified'
            err_note = """'You must set a 'Pipeline Scripts Directory' to run the\
             pipeline Perl scripts. Please go into the admin control panel and add\
             a 'Pipeline Scripts Directory'.  Talk to your administrator for details."""
            self.record_error(err_title, err_note)
            return None

        return script_directory_info

    def run_script_on_dataset(self):
        """
        Using dataset information, decide whether to run:
            (1)
            (2)
        """

        # (1) Get the script directory
        #
        script_directory_info = self.get_script_directory_info()
        if script_directory_info is None:
            return (False, self.get_err_msg_obj())


        if self.dataset.is_vcf_file():     # (2a) bsub command for a VCF file
            return self.run_vcf_file_script(script_directory_info)

    def run_vcf_file_script(self, script_directory_info):
        """
        Run the bsub command for a VCF file.
        e.g. bsub perl analyseVCF.pl (directory name)
        """
        if self.err_found:
            return (False, self.get_err_msg_obj())

        # (1) Make sure the 'analyseVCF.pl' command is in the
        #   specified 'Pipeline Scripts Directory'
        #
        script_cmd = join(script_directory_info.script_directory, 'analyseVCF.pl' )
        if not isfile(script_cmd):
            err_title='The "analyseVCF.pl" file was not found'
            err_note = """The "analyseVCF.pl" file was not found at this location: \n%s\n
            To change the path to the script, please go into the admin control\
             panel and modify the 'Pipeline Scripts Directory'.\n\
            Talk to your administrator for details.""" % (script_cmd)
            err_msg_obj = self.record_error(err_title, err_note)
            return (False, err_msg_obj)

        # Format the full bsub command with target containing
        #   input files
        #
        command_str = 'bsub perl {0} {1}'.format(script_cmd,\
                dataset.file_directory)


        # (2) Create a run object and save the command being run
        #
        dsr = DatasetScriptRun(dataset=dataset)
        dsr.notes = command_str
        dsr.save()

        # (3) Update the dataset status to 'in process'
        #
        dataset.set_status_processing_started()

        # (4) Run the script -- not waiting for output
        #
        full_args = command_str.split()
        print ('full_args', full_args)
        run_script(full_args)

        return (True, dsr)
