from __future__ import print_function
import os, sys
from os.path import dirname, realpath

if __name__=='__main__':
    django_dir = dirname(dirname(dirname(realpath(__file__))))
    sys.path.append(django_dir)
    os.environ['DJANGO_SETTINGS_MODULE'] = 'tb_website.settings.local'

    # Allows the working environ to get set-up, apps registered, etc
    #
    import django
    django.setup()

from datetime import datetime
import requests

from apps.predict.models import PredictDataset,\
            DATASET_STATUS_FILE_RETRIEVAL_COMPLETE
from apps.predict.script_run_helper import run_script_on_dataset


class PipelineScriptRunner:


    @staticmethod
    def run_datasets():
        """
        Run datasets through the pipeline whose files have just been retrieved.

        Note: These run serially--one after another--if you have several datasets, run them one at a time.

        Try the "PipelineScriptRunner.run_next_dataset()"
        """
        qset = PredictDataset.objects.filter(status=DATASET_STATUS_FILE_RETRIEVAL_COMPLETE)

        print ('Found {0} dataset(s) to run through pipeline'.format(qset.count()))

        for pd in qset:
            run_script_on_dataset(pd)

    @staticmethod
    def run_next_dataset():
        """
        Run datasets through the pipeline whose files have just been retrieved.

        Note: These run serially--one after another--if you have several datasets, run them one at a time.

        Try the "PipelineScriptRunner.run_next_dataset()"
        """
        pd = PredictDataset.objects.filter(status=DATASET_STATUS_FILE_RETRIEVAL_COMPLETE).first()
        if pd is None:
            print('No dataset found to run through pipeline')
            return

        run_script_on_dataset(pd)

    @staticmethod
    def run_specific_dataset(dataset_id):
        """
        Run datasets through the pipeline whose files have just been retrieved.

        Note: These run serially--one after another--if you have several datasets, run them one at a time.

        Try the "PipelineScriptRunner.run_next_dataset()"
        """
        # Get PredictDataset
        #
        try:
            dataset = PredictDataset.objects.get(pk=dataset_id)
        except PredictDataset.DoesNotExist:
            print ('Failed.  There is not "PredictDataset" with db id: {0}'.format(dataset_id))
            return False

        if dataset.status < DATASET_STATUS_FILE_RETRIEVAL_COMPLETE:
            print ('Failed.  Files not yet available for this "PredictDataset" db id: {0}'.format(dataset_id))
            return False

        run_script_on_dataset(dataset)



if __name__=='__main__':
    args = sys.argv
    if len(args) == 1:

        # PredictDatasets where files have been downloaded
        #  and processing has NOT been attempted
        #
        PipelineScriptRunner.run_datasets()

    elif len(args) == 2 and args[1] == '--next':

        # The next PredictDataset where files have been downloaded
        #  and processing has NOT been attempted
        #
        PipelineScriptRunner.run_next_dataset()

    elif len(args) == 2 and args[1].isdigit():

        # Single PredictDataset specified by database id
        #   regardless of processing status
        #
        PipelineScriptRunner.run_specific_dataset(dataset_id=args[1])

    else:
        print('-' * 40)
        print("""Regular run of new dropbox links:
    >python pipeline_script_runner.py

Retry dropbox links with errors:
    >python pipeline_script_runner.py --next

Retrieve dropbox link files for a specifiec PredictDataset:
    >python pipeline_script_runner.py (dataset id)
   e.g. >python pipeline_script_runner.py 102
        """)
