from apps.predict.models import DropboxRetrievalLog#PredictDatasetStatus, PredictDataset, PredictDatasetNote, DatasetScriptRun, ScriptToRun
from apps.predict.dropbox_retriever import DropboxRetriever
from apps.script_helper.script_runner_basic import run_script

class ErrMsg:
    def __init__(self, title, note):
        self.title = title
        self.note = note

def run_dropbox_retrieval(dbox_log):
    assert isinstance(dbox_log_obj, DropboxRetrievalLog), 'dbox_log must be a DropboxRetrievalLog object'

    """
    dr = DropboxRetriever(example_dlink, dest_dir, GENTB_FILE_PATTERNS)
    if dr.err_found:
        print dr.err_msg
        sys.exit(1)

    # Get the metadata
    #
    if not dr.step1_retrieve_metadata():
        print dr.err_msg
        sys.exit(1)

    # Does it have what we want?
    #
    if not dr.step2_check_file_matches():
        print dr.err_msg
        sys.exit(1)

    print dr.matching_files_metadata
    #sys.exit(1)

    # Download the files
    #
    if not dr.step3_retrieve_files():
        print dr.err_msg
        sys.exit(1)
"""


    # (2) Create a run object
    #
    dsr = DatasetScriptRun(dataset=dataset)
    dsr.save()

    # (2) Format actual command line args
    #
    json_args = dataset.get_script_args_json(run_md5=dsr.md5, as_list=True)

    # (3) Combine script + dataset args
    #
    full_args = chosen_script.get_script_args_as_list() + json_args
    print ('full_args', full_args)

    # (4) Save the command being run
    #
    dsr.notes='%s' % ' '.join(full_args)
    dsr.save()

    # (5) Update the dataset status to 'in process'
    #
    dataset.set_status_processing_started()

    # (6) Run the script -- not waiting for output
    #
    run_script(full_args)


    return (True, dsr)
