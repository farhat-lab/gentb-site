import sys
from apps.dropbox_helper.models import DropboxRetrievalLog
from apps.predict.dropbox_retriever import DropboxRetriever
from apps.script_helper.script_runner_basic import run_script
from apps.dropbox_helper.forms import DropboxRetrievalParamsForm

class ErrMsg:
    def __init__(self, title, note):
        self.title = title
        self.note = note

class DropboxRetrievalRunner:

    def __init__(self, json_args_string):
        assert json_args_string is not None and len(json_args_string) > 0,\
            "json_args_string must contain arguments that make the DropboxRetrievalParamsForm valid"

        # Error messages
        self.err_found = False
        self.err_msg = None

        # Basic params
        self.dropbox_url = None
        self.destination_directory = None
        self.callback_url = None
        self.callback_md5 = None

        # Set params
        if not self.set_params(json_args_string):
            return

        # Run retrieval_error
        self.run_dropbox_retrieval()

    def set_error_message(self, m):
        self.err_found = True
        self.err_msg = m

    def set_params(self, json_args_string)
        try:
            json_args = json.loads(json_args_string)
        except:
            # MAJOR ERROR!  LOG THIS!  CANNOT MAKE CALLBACK
            self.set_error_message("These arguments could NOT be converted to JSON: {0}".format(json_args_string))
            return False

        # Validate the params
        #
        f = DropboxRetrievalParamsForm(json_args)
        if not f.is_valid():
            # MAJOR ERROR!  LOG THIS!  LIKELY CANNOT MAKE CALLBACK
            self.set_error_message("Invalid arguments: {0}".format(f.errors))
            return False

        self.dropbox_url = f.cleaned_data['dropbox_url']
        self.destination_directory = f.cleaned_data['destination_directory']
        self.callback_url = f.cleaned_data['callback_url']
        self.callback_md5 = f.cleaned_data['callback_md5']


    def send_callback_message(self, success_flag, result_data):

        params = dict(success=success_flag,
                    callback_md5=self.callback_md5,
                    result_data=result_data
                    )

        # Make the callback update
        #
        r = requests.post(self.callback_url,
                        data=params,
                        )

        print r.status_code
        print r.text
        if not r.status_code = 200:
            # LOG THIS
            print 'FAIL'


    def run_dropbox_retrieval(json_args_string):
        if self.err_found:
            return False

        dr = DropboxRetriever(self.dropbox_url,
                        self.destination_directory)

        if dr.err_found:
            send_callback_message(False, dr.get_err_msg_as_dict())
            return False

        # Get the metadata
        #
        if not dr.step1_retrieve_metadata():
            send_callback_message(False, dr.get_err_msg_as_dict())
            return False

        # Does it have what we want?
        #
        if not dr.step2_check_file_matches():
            send_callback_message(False, dr.get_err_msg_as_dict())
            return False

        # Download the files
        #
        if not dr.step3_retrieve_files():
            send_callback_message(False, dr.get_err_msg_as_dict())
            return False

        send_callback_message(True, dr.final_file_paths)


if __name__=='__main__':
    args = sys.argv
    if len(args) != 2:
        print """>python drobox_retrieval_runner.py '{ json string }'"""
        print """Example:\npython drobox_retrieval_runner.py {"dropbox_url": "https://a-dropbox-shared-link.com", "callback_md5": "8fec9fefa93095fc94a68f495e24325b", "destination_directory": "/an-existing-dir-to-put-files", "callback_url": "https://myserver.com/predict/file-retrieval-results"}"""
    else:
        dr = DropboxRetrievalRunner(sys.argv[1])
