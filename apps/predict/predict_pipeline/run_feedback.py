"""
This script checks the local "output" directory and sends feedback to the server
about the job.

Dataset file diretory/
        - feedback.conf
        - /output/
            - result.json
            - matrix.csv
"""
from __future__ import print_function
from os.path import dirname, join, isdir, isfile, getsize, realpath
import os, sys
import urllib

"""
# example CALLBACK_INFO_DICT
CALLBACK_INFO_DICT = {"file_directory": "/home/gentb_test/tbdata_00000112",
"run_md5": "bb897a28f59f93ad115a6faa42f4918d",
"admin_url": "https://gentb.hms.harvard.edu/gentb-admin/predict/predictdataset/1/",
"callback_url": "http://127.0.0.1:8000/predict/pipeline-run-results-notice/",
#"callback_url": "https://gentb.hms.harvard.edu/predict/pipeline-run-results-notice/",
"dataset_id": 112,
"user_email": "tbuser@harvard.edu"}
"""
with open(join(sys.argv[1], 'feedback.conf'), 'r') as fhl:
    FEEDBACK_CONF = json.loads(fhl.read())

CALLBACK_INFO_DICT = FEEBACK_CONF['callback_info_dict']


class GenTBStatusFeedback:

    CURRENT_DIRECTORY = dirname(realpath(__file__))
    OUTPUT_DIRECTORY = join(CURRENT_DIRECTORY, FEEDBACK_CONF['RESULT_OUTPUT_DIRECTORY_NAME'])

    # note file names originate in utils.result_file_info.py
    #
    RESULTS_FILE_FULLPATH = join(OUTPUT_DIRECTORY, FEEDBACK_CONF['RESULT_JSON_FILE_NAME'])
    MATRIX_FILE_FULLPATH = join(OUTPUT_DIRECTORY, FEEDBACK_CONF['MATRIX_CSV_FILE_NAME'])
    HEATMAP_FILE_FULLPATH = join(OUTPUT_DIRECTORY, FEEDBACK_CONF['HEATMAP_HTML_FILE_NAME'])

    FILE_FULLPATH_NAMES = [RESULTS_FILE_FULLPATH,\
                MATRIX_FILE_FULLPATH,\
                HEATMAP_FILE_FULLPATH]

    # Ex: ('results file', 'matrix file', 'heatmap file')
    #
    EXPECTED_FILE_NAME_LIST = FEEDBACK_CONF['EXPECTED_FILE_DESCRIPTIONS']

    EXPECTED_FILE_DESCRIPTIONS = zip(FILE_FULLPATH_NAMES, EXPECTED_FILE_NAME_LIST)

    def __init__(self):
        self.job_checked = False
        self.success = False
        self.error_message = None

        self.run_job_result_check()
        print ('did_job_succceed', self.did_job_succceed())
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
            "Do not call this before running 'run_job_result_check()'"

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
        fh = open(self.RESULTS_FILE_FULLPATH, 'r')
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

        print ('Callback url: {0}\nParams: {1}'.format(\
                    callback_url, callback_dict))

        callback_params = urllib.urlencode(callback_dict)
        f = urllib.urlopen(callback_url, callback_params)
        print ('Callback result: {0}'.format(f.read()))



    def run_job_result_check(self):
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
        # Check each result file
        #   - Does it exist?
        #   - Does it have data?
        # ------------------------
        for fullname, human_name in self.EXPECTED_FILE_DESCRIPTIONS:

            if not isfile(fullname):
                self.add_error("The %s was not found: %s"\
                        % (human_name, fullname))
                return

            if getsize(fullname) == 0:
                self.add_error("The %s was empty: %s"\
                        %  (human_name, fullname))
                return

        self.success = True

if __name__ == '__main__':
    gentb_status = GenTBStatusFeedback()
