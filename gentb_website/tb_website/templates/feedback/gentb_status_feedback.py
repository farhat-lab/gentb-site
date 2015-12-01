"""
This script checks the local "output" directory and sends feedback to the server
about the job.

/output/
    - result.json
    - matrix.csv
"""
from os.path import dirname, join, isdir, isfile, getsize, realpath
import os, sys
import urllib


# example CALLBACK_INFO_DICT
CALLBACK_INFO_DICT = {"file_directory": "/home/gentb_test/tbdata_00000112",
"run_md5": "bb897a28f59f93ad115a6faa42f4918d",
"admin_url": "https://gentb.hms.harvard.edu/gentb-admin/predict/predictdataset/1/",
"callback_url": "http://127.0.0.1:8000/predict/pipeline-run-results-notice/",
#"callback_url": "https://gentb.hms.harvard.edu/predict/pipeline-run-results-notice/",
"dataset_id": 112,
"user_email": "tbuser@harvard.edu"}

#CALLBACK_INFO_DICT = {{ callback_info_dict }}


class GenTBStatusFeedback:

    CURRENT_DIRECTORY = dirname(realpath(__file__))
    OUTPUT_DIRECTORY = join(CURRENT_DIRECTORY, 'output')
    RESULTS_FILE_NAME = join(OUTPUT_DIRECTORY, 'result.json')
    MATRIX_FILE_NAME = join(OUTPUT_DIRECTORY, 'matrix.csv')

    def __init__(self):
        self.job_checked = False
        self.success = False
        self.error_message = None

        self.check_if_job_succeeded()
        print 'did_job_succceed', self.did_job_succceed()
        self.send_feedback_to_gentb()

    def add_error(self, msg):
        """
        Store error messages
        """
        self.success = False
        self.error_message = msg

    def did_job_succceed(self):
        """
        Was the expected output found?
        """
        assert self.job_checked is True,\
            "Do not call this before running 'check_if_job_succeeded()'"

        return self.success


    def get_callback_params(self):
        """
        Format the callback parameters
            run_md5 - identify job
            success - did it work
            result_data -   Success: data in results.json
                            Fail: error message

        """
        callback_params = dict(run_md5=CALLBACK_INFO_DICT.get('run_md5', None))

        # --------------------
        # Job Failed
        # --------------------
        if not self.did_job_succceed():
            callback_params['success'] = False
            callback_params['result_data'] = self.error_message
            return callback_params

        # --------------------
        # Job succeeded
        # --------------------
        callback_params['success'] = True

        # Send back results.json data
        #
        fh = open(self.RESULTS_FILE_NAME, 'r')
        result_info = fh.read()
        fh.close()
        callback_params['result_data'] = result_info

        return callback_params


    def send_feedback_to_gentb(self):
        """
        Using the callback_url in CALLBACK_INFO_DICT, notify the
        genTB web server if the job ran successfully
        """
        callback_url = CALLBACK_INFO_DICT.get('callback_url', None)
        if callback_url is None:
            return

        callback_dict = self.get_callback_params()

        callback_params = urllib.urlencode(callback_dict)
        f = urllib.urlopen(callback_url, callback_params)
        print f.read()



    def check_if_job_succeeded(self):
        """
        Check if the expected directory/file output is found
        """
        if self.job_checked is True:
            return

        self.job_checked = True

        # Does the output directory exist:
        #
        if not isdir(self.OUTPUT_DIRECTORY):
            self.add_error("The output directory does not exist: %s"\
                % self.OUTPUT_DIRECTORY)
            return

        # ------------------------
        # Check results file
        #   - Does it exist?
        #   - Does it have data?
        # ------------------------
        if not isfile(self.RESULTS_FILE_NAME):
            self.add_error("The results file was not found: %s"\
                % self.RESULTS_FILE_NAME)
            return

        if getsize(self.RESULTS_FILE_NAME) == 0:
            self.add_error("The results file was empty: %s"\
                % self.RESULTS_FILE_NAME)
            return

        # ------------------------
        # Check matrix file
        #   - Does it exist?
        #   - Does it have data?
        # ------------------------
        if not isfile(self.MATRIX_FILE_NAME):
            self.add_error("The matrix file was not found: %s"\
                % self.MATRIX_FILE_NAME)
            return

        if getsize(self.MATRIX_FILE_NAME) == 0:
            self.add_error("The matrix file was empty: %s"\
                % self.MATRIX_FILE_NAME)
            return

        self.success = True

if __name__ == '__main__':
    gentb_status = GenTBStatusFeedback()
