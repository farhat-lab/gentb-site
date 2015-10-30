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

from apps.dropbox_helper.models import DropboxRetrievalLog
from apps.dropbox_helper.dropbox_retriever import DropboxRetriever
from apps.script_helper.script_runner_basic import run_script
from apps.dropbox_helper.forms import DropboxRetrievalParamsForm


class DropboxRetrievalRunner:

    def __init__(self, dbox_log):
        assert isinstance(dbox_log, DropboxRetrievalLog), "dbox_log must be an instance of DropboxRetrievalLog"

        self.dbox_log = dbox_log
        self.predict_dataset = self.dbox_log.dataset

        # Error messages
        self.err_found = False
        self.err_msg = None

        # Run retrieval_error
        self.run_dropbox_retrieval()

    def set_error_message(self, m):
        self.err_found = True
        self.err_msg = m


    def record_retrieval_error(self, error_msg):

        self.set_error_message(error_msg)
        #msg = "There was an error retrieving the dropbox files."

        # Update dataset status
        #
        self.predict_dataset.set_status_file_retrieval_error()

        # Update dropbox log
        #
        self.dbox_log.retrieval_error = error_msg
        self.dbox_log.set_retrieval_end_time(self, files_retrieved=False)
        self.dbox_log.save()

    def run_dropbox_retrieval(self):

        if self.err_found:
            return False


        # Reset log attributes and mark the retrieval start time
        #
        self.dbox_log.set_retrieval_start_time(with_reset=True)
        self.dbox_log.save()

        # Update PredictDataset status to 'file retrieval started'
        #
        self.dbox_log.dataset.set_status_file_retrieval_started()

        dr = DropboxRetriever(self.predict_dataset.dropbox_url,
                        self.predict_dataset.file_directory)

        if dr.err_found:
            self.record_retrieval_error(dr.err_msg)
            return False

        # Get the metadata
        #
        if not dr.step1_retrieve_metadata():
            self.record_retrieval_error(dr.err_msg)
            return False

        # Does it have what we want?
        #
        if not dr.step2_check_file_matches():
            self.record_retrieval_error(dr.err_msg)
            return False

        # Download the files
        #
        if not dr.step3_retrieve_files():
            self.record_retrieval_error(dr.err_msg)
            return False

        # ----------------------
        # Success!
        # ----------------------
        self.predict_dataset.set_status_file_retrieval_complete()
        self.dbox_log.selected_files = dr.final_file_paths
        self.dbox_log.set_retrieval_end_time(files_retrieved=True)
        self.dbox_log.save()

    @staticmethod
    def retrieve_new_dropbox_files(**kwargs):
        """
        Are there DropboxRetrievalLog objects
        where the download process hasn't started?

        Yes, start retrieving the files

        Alternate: **kwargs={ 'retry_files_with_errors' : True}
            Try to download files where the process encountered an error
        """

        retry_files_with_errors = kwargs.get('retry_files_with_errors', True)
        if retry_files_with_errors is True:
            # Files where download encountered an error
            #
            qs_args = dict(retrieval_error__isnull=False)
        else:
            # Files where download process hasn't started
            #
            qs_args = dict(retrieval_start__isnull=True)

        # Get files that haven't been retrieved from dropbox urls
        #
        dbox_logs_to_check = DropboxRetrievalLog.objects.filter(**qs_args\
                            ).exclude(files_retrieved=True)

        cnt = dbox_logs_to_check.count()
        if cnt == 0:
            print 'All set.  Nothing to check'
            return

        print 'Checking {0} link(s)'.format(cnt)
        
        for dbox_log in dbox_logs_to_check:

            # Go get the files!
            print 'run dlog', dbox_log
            print 'with params:', dbox_log.get_dropbox_retrieval_script_params()
            dr = DropboxRetrievalRunner(dbox_log)


if __name__=='__main__':
    args = sys.argv
    if len(args) == 1:
        DropboxRetrievalRunner.retrieve_new_dropbox_files()

    elif len(args) == 2 and args[1] == '--retry':
        retry_param= dict(retry_files_with_errors=True)
        DropboxRetrievalRunner.retrieve_new_dropbox_files(**retry_param)

    else:
        print '-' * 40
        print """Regular run of new dropbox links:
    >python dropbox_retrieval_runner.py

Retry dropbox links with errors:
    >python dropbox_retrieval_runner.py --retry
        """
