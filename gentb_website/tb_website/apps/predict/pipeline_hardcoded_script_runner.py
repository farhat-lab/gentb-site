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

VCF_ANALYSIS_SCRIPT = 'analyseVCF.pl'
FASTQ_ANALYSIS_SCRIPT = 'analyseNGS.pl'


class PipelineScriptRunner(object):

    def __init__(self, dataset):
        assert isinstance(dataset, PredictDataset),\
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

    def run_script_on_dataset(self):

        script_directory = self.step1_get_script_directory_info()
        if script_directory is None:
            return (False, self.get_err_msg_obj())

        script_command = self.step2_get_script_command(script_directory)
        if script_command is None:
            return (False, self.get_err_msg_obj())

        return self.step3_run_command(script_command)


    def step1_get_script_directory_info(self):
        """
        Retrieve PipelineScriptsDirectory object from the database
        """
        script_directory_info = PipelineScriptsDirectory.objects.filter(is_chosen_directory=True).first()
        if script_directory_info is None:
            err_title='No pipeline directory specified'
            err_note = """'You must set a 'Pipeline Scripts Directory' to run the\
             pipeline Perl scripts. Please go into the admin control panel and add\
             a 'Pipeline Scripts Directory'.  Talk to your administrator for details."""
            self.record_error(err_title, err_note)
            return None

        return script_directory_info

    def step2_get_script_command(self, script_directory_info):
        """
        Using dataset information, to decide whether to run:
            (1) script for a VCF file
            (2) script for FastQ files
        """
        if script_directory_info is None:
            return None

        # Formate either a VCF or FastQ pipeline command
        #
        if self.dataset.is_vcf_file():     # (2a) bsub command for a VCF file
            command_to_run = self.get_vcf_script_command(script_directory_info)

        elif self.dataset.is_fastq_file(): # (2b) bsub command for FastQ files
            command_to_run = self.get_fastq_script_command(script_directory_info)

        else:
            err_title='Not VCF or FastQ file'
            err_note = 'Could not determine the file type.\
             Database contained: "%s"' % (dataset.file_type)
            err_msg_obj = self.record_error(err_title, err_note)
            return None

        return command_to_run


    def step3_run_command(self, command_to_run):

        if command_to_run is None:
            return (False, self.get_err_msg_obj())

        # (1) Create a run object and save the command being run
        #
        dsr = DatasetScriptRun(dataset=self.dataset)
        dsr.notes = command_to_run
        dsr.save()

        # (2) Update the dataset status to 'in process'
        #
        self.dataset.set_status_processing_started()

        # (3) Run the script -- not waiting for output
        #
        full_args = command_to_run.split()
        print ('full_args', full_args)
        run_script(full_args)

        return (True, dsr)


    def get_fastq_script_command(self, script_directory_info):
        """
        Run the bsub command for FastQ files.
        e.g. bsub perl analyseNGS.pl (directory name)
        """
        if self.err_found:
            return None

        # (1) Make sure the 'analyseVCF.pl' command is in the
        #   specified 'Pipeline Scripts Directory'
        #
        script_cmd = join(script_directory_info.script_directory, FASTQ_ANALYSIS_SCRIPT)
        if not isfile(script_cmd):
            err_title='The "%s" file was not found' % (FASTQ_ANALYSIS_SCRIPT)
            err_note = """The "%s" file was not found at this location: \n%s\n
            To change the path to the script, please go into the admin control\
             panel and modify the 'Pipeline Scripts Directory'.\n\
            Talk to your administrator for details.""" % (FASTQ_ANALYSIS_SCRIPT, script_cmd)
            err_msg_obj = self.record_error(err_title, err_note)
            return None

        # Format the full bsub command with target containing
        #   input files
        #
        if self.dataset.is_fastq_single_ended():
            command_str = 'bsub perl {0} 0 . {1}'.format(script_cmd,\
                dataset.file_directory)
        elif self.dataset.is_fastq_pair_ended():
            command_str = 'bsub perl {0} 1 . {1}'.format(script_cmd,\
                dataset.file_directory)
        else:
            err_title = 'FastQ: single-ended or pair-ended?'
            err_note = 'Could not determine single-ended or pair-ended FastQ type.\
             Database contained: "%s"' % (dataset.fastq_type)
            err_msg_obj = self.record_error(err_title, err_note)
            return None

        return command_str



    def get_vcf_script_command(self, script_directory_info):
        """
        Run the bsub command for a VCF file.
        e.g. bsub perl analyseVCF.pl (directory name)
        """
        if self.err_found:
            return None

        # (1) Make sure the 'analyseVCF.pl' command is in the
        #   specified 'Pipeline Scripts Directory'
        #
        script_cmd = join(script_directory_info.script_directory, VCF_ANALYSIS_SCRIPT )
        if not isfile(script_cmd):
            err_title='The "%s" file was not found' % (VCF_ANALYSIS_SCRIPT)
            err_note = """The "%s" file was not found at this location: \n%s\n
            To change the path to the script, please go into the admin control\
             panel and modify the 'Pipeline Scripts Directory'.\n\
            Talk to your administrator for details.""" % (VCF_ANALYSIS_SCRIPT, script_cmd)
            err_msg_obj = self.record_error(err_title, err_note)
            return None

        # Format the full bsub command with target containing
        #   input files
        #
        command_str = 'bsub perl {0} {1}'.format(script_cmd,\
                self.dataset.file_directory)

        return command_str

if __name__ == '__main__':
    pass
    """
    dataset = PredictDataset.objects.first()

    # get some Dataset
    pipeline_runner = PipelineScriptRunner(dataset)


    pipeline_runner.run_script_on_dataset()

    # OR

    script_directory = pipeline_runner.step1_get_script_directory_info()
    if script_directory is not None:
        script_command = pipeline_runner.step2_get_script_command(script_directory)
        if script_command is not NOne:
            pipeline_runner.step3_run_command(script_command)


    """
