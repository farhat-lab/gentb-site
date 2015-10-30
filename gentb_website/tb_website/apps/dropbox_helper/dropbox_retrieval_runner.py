import sys
from apps.dropbox_helper.models import DropboxRetrievalLog
from apps.predict.dropbox_retriever import DropboxRetriever
from apps.script_helper.script_runner_basic import run_script
from apps.dropbox_helper.forms import DropboxRetrievalParamsForm

class ErrMsg:
    def __init__(self, title, note):
        self.title = title
        self.note = note

def run_dropbox_retrieval(**kwargs):
    #assert isinstance(dbox_log_obj, DropboxRetrievalLog),\
    #    'dbox_log must be a DropboxRetrievalLog object'
    assert kwargs is not None and len(kwargs) > 0,\
        "kwargs must contain arguments that make the DropboxRetrievalParamsForm valid"

    # Validate the params
    #
    f = DropboxRetrievalParamsForm(kwargs)
    if not f.is_valid():
        # MAJOR ERROR!  LOG THIS!
        raise LookupError("Invalid form: {0}".format(f.errors))

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

if __name__=='__main__':
    args = sys.argv
    if len(args) != 2:
        print """>python drobox_retrieval_runner.py '{ json string }'"""
        print """Example:\npython drobox_retrieval_runner.py """
