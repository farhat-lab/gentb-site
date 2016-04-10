"""
Examine existing PredictDataset objects.

- Retrieve files from confirmed dropbox links

"""

from django.core.management.base import BaseCommand, CommandError

from apps.predict.models import PredictDataset
from apps.dropbox.models import DropboxRetrievalLog
from apps.dropbox.retriever import DropboxRetriever

import logging
LOGGER = logging.getLogger('apps.dropbox.management.commands.dropbox_retrieval')

class Command(BaseCommand):
    help = """Regular run of new dropbox links:

  manage.py dropbox_retrieval

Retry dropbox links with errors:

  manage.py dropbox_retrieval --retry

Retrieve dropbox link files for a specifiec PredictDataset:

  manage.py dropbox_retrieval <dataset_id>

"""

    def add_arguments(self, parser):
        parser.add_argument('--dataset_id', nargs='+', type=int,
            help='Retrieve dropbox link files for a specifiec PredictDataset')

        parser.add_argument('--retry',
            action='store_true',
            dest='retry',
            default=False,
            help='Retry dropbox links with errors')

    def handle(self, dataset_id=None, retry=False, **options):
        if dataset_id is not None:
            # Single PredictDataset specified by database id, regardless of status
            self.retrieve_dataset_files(*dataset_id)
        elif retry:
            # PredictDatasets where download previously failed
            self.retrieve_new_dropbox_files(with_errors=True)
        else:
            self.retrieve_new_dropbox_files()

    def retrieve_dataset_files(self, *ids):
        """
        Given a dataset id, retrieve the dropbox link files, REGARDLESS of previous attempts

        For running form the command line, output is print statments
        """
        # Get PredictDataset and DropboxRetrievalLog
        #
        pds = PredictDataset.objects.filter(pk__in=ids)
        logs = DropboxRetrievalLog.objects.filter(dataset_id__in=ids)

        for pks in (pds, logs):
            not_pks = set(ids) ^ set(pks.values_list('id', flat=True))
            if not_pks:
                name = pks.model.__name__
                items = ", ".join([str(pk) for pk in not_pks])
                raise CommandError("No '%s' with ids %s" % (name, items))

        # Get first matching DropboxRetrievalLog object for retrieval
        #
        for pk in pds:
            log = DropboxRetrievalLog.objects.filter(dataset_id=pk).first()
            dr = DropboxRetrievalRunner(log)

    def retrieve_new_dropbox_files(self, with_errors=True):
        """
        Are there DropboxRetrievalLog objects
        where the download process hasn't started?

        Yes, start retrieving the files

        Alternate: **kwargs={ 'retry_files_with_errors' : True}
            Try to download files where the process encountered an error
        """
        LOGGER.debug("retrieve_new_dropbox_files")

        qs = DropboxRetrievalLog.objects.exclude(files_retrieved=True)
        if with_errors is True:
            # Files where download encountered an error
            #
            qs = qs.filter(retrieval_error__isnull=False)
        else:
            # Files where download process hasn't started
            #
            qs = qs.filter(retrieval_start__isnull=True)

        if qs.count() == 0:
            status_msg = 'All set.  Nothing to check'
            print(status_msg)
            LOGGER.debug(status_msg)
            return

        print('Checking {0} link(s)'.format(qs.count()))

        for dbox_log in qs:
            print('Get files for: ', dbox_log)
            dr = DropboxRetrievalRunner(dbox_log)



class DropboxRetrievalRunner(object):

    def __init__(self, dbox_log):
        assert isinstance(dbox_log, DropboxRetrievalLog), "dbox_log must be an instance of DropboxRetrievalLog"

        LOGGER.debug("Initiate dropbox retrieval for: %s", dbox_log.dataset)
        LOGGER.debug(" - dropbox_url:  dropbox retrieval for: %s",\
            dbox_log.dataset.dropbox_url)

        self.dbox_log = dbox_log
        self.predict_dataset = self.dbox_log.dataset

        # Error messages
        self.err_found = False
        self.err_msg = None

        # Run retrieval_error
        self.run_dropbox_retrieval()

    def set_error_message(self, m):
        """
        Set internal error flag and store the message within this object
        """
        LOGGER.error("Error found during dropbox file retrieval: %s", m)

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
                        self.predict_dataset.file_directory,
                        file_patterns=self.predict_dataset.get_file_patterns())

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

