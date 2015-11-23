from os.path import join, isdir
import os
from apps.dropbox_helper.dropbox_retriever import DropboxRetriever
from apps.dropbox_helper.models import DropboxRetrievalLog
from django.conf import settings


def get_dropbox_metadata(predict_dataset):
    """
    Wrap the DropboxRetriever function
    - (True, DropboxRetrievalLog object)
    - (False, error message string)
    """
    if predict_dataset is None:
        return (False, "The dataset was not found.")


    # Initialize
    #
    dr = DropboxRetriever(predict_dataset.dropbox_url,
                          destination_dir=predict_dataset.file_directory,
                          file_patterns=self.predict_dataset.get_file_patterns())


    db_log = DropboxRetrievalLog(dataset=predict_dataset)

    if dr.err_found:
        db_log.file_metadata_err_msg = dr.err_msg
        db_log.save()
        return (False, dr.err_msg)

    # Get the metadata
    #
    if not dr.step1_retrieve_metadata():
        db_log.file_metadata_err_msg = dr.err_msg
        db_log.save()
        return (False, dr.err_msg)

    # Does it have what we want?
    #
    if not dr.step2_check_file_matches():
        db_log.file_metadata_err_msg = dr.err_msg
        db_log.save()
        return (False, dr.err_msg)

    # Yes!
    db_log.file_metadata = dr.matching_files_metadata
    db_log.save()
    return (True, dr)


def get_dropbox_metadata_from_link(dropbox_link, file_patterns=None):
    """
    Wrap the DropboxRetriever function
    - (True, metadata-selected files)
    - (False, error message string)
    """
    if dropbox_link is None:
        return (False, "The dataset was not found.")

    # This directory doesn't actually get used
    #
    tmp_dir = join(settings.TB_SHARED_DATAFILE_DIRECTORY, 'tmp')
    if not isdir(tmp_dir):
        os.makedirs(tmp_dir)

    # Initialize
    #
    if file_patterns:
        dr = DropboxRetriever(dropbox_link,\
                  destination_dir=tmp_dir,\
                  file_patterns=file_patterns)

    else:
        dr = DropboxRetriever(dropbox_link,\
                          destination_dir=tmp_dir)

    if dr.err_found:
        return (False, dr.err_msg)

    # Get the metadata
    #
    if not dr.step1_retrieve_metadata():
        return (False, dr.err_msg)

    # Does it have what we want?
    #
    if not dr.step2_check_file_matches():
        return (False, dr.err_msg)

    # Yes!
    return (True, dr.matching_files_metadata)
