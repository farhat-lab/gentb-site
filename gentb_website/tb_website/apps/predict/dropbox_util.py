from apps.predict.dropbox_retriever import DropboxRetriever
from apps.predict.models import DropboxRetrievalLog



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
                          destination_dir=predict_dataset.file_directory)

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