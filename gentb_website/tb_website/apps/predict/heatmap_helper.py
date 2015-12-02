"""
Given a dataset, retrieve and format its heatmap HTML
for use in a template
"""
from os.path import isdir, isfile, join, getsize
from apps.utils.result_file_info import HEATMAP_HTML_FILE_NAME,\
        RESULT_OUTPUT_DIRECTORY_NAME
from apps.predict.models import PredictDataset

import logging
LOGGER = logging.getLogger(__name__)

class HeatmapHelper(object):

    LINES_TO_REMOVE = """<!DOCTYPE html>
    <html>
    <head>
    <meta charset="utf-8" />
    </head>
    <body style="background-color:white;">
    </body>
    </html>""".split('\n')
    LINES_TO_REMOVE = [x.strip() for x in LINES_TO_REMOVE]

    def __init__(self, dataset):
        """
        For a given PredictDataset, retrieve the related heatmap html
        """
        self.dataset = dataset
        self.heatmap_html_lines = None
        self.heatmap_html = None

        self.has_error = False
        self.error_message = None

        self.check_dataset()
        if not self.has_error:
            self.load_heatmap_html()

    def add_error(self, msg):
        """
        Store error messages
        """
        self.has_error = True
        self.error_message = msg

    def check_dataset(self):
        if not isinstance(self.dataset, PredictDataset):
            self.add_error("The dataset must be a PredictDataset object")
            return

        if not self.dataset.has_prediction:
            self.add_error("Sorry, this dataset does NOT have a prediction")
            return

    def load_heatmap_html(self):
        if self.has_error:
            return

        # get heatmap file name
        #
        heatmap_file_fullname = join(self.dataset.file_directory,\
                            RESULT_OUTPUT_DIRECTORY_NAME,\
                            HEATMAP_HTML_FILE_NAME)

        # Is it a file?
        #
        if not isfile(heatmap_file_fullname):
            LOGGER.error("The heatmap file for dataset id %d was not found: %s"\
                        % (self.dataset.id, heatmap_file_fullname))
            self.add_error("Sorry, the heatmap file was not found.")
            return

        # Is the size greater than zero?
        #
        if getsize(heatmap_file_fullname) == 0:
            LOGGER.error("Empty heatmap file for dataset id %d: %s"\
            % (self.dataset.id, heatmap_file_fullname))
            err_msg = "Sorry, the heatmap file was empty."
            return

        # Try to open it and read the HTML
        #
        try:
            fh = open(heatmap_file_fullname, 'r')
            self.heatmap_html_lines = fh.readlines()
            fh.close()
        except Exception as e:
            LOGGER.error("Error trying to open the heatmap file for dataset id %d: %s\n\n%s"\
                        % (self.dataset.id, heatmap_file_fullname, e))
            err_msg = "Sorry, an error was found when opening the heatmap file."
            return

        if self.heatmap_html_lines is None or len(self.heatmap_html_lines) == 0:
            LOGGER.error("Error no content was read from the heatmap file for dataset id %d: %s\n\n%s"\
                        % (self.dataset.id, heatmap_file_fullname, e))
            err_msg = "Sorry, no content was read from the heatmap file."

        self.format_heatmap()

    def is_line_to_remove(self, heatmap_line):
        if heatmap_line is None:
            return None

        for line_to_go in self.LINES_TO_REMOVE:
            print 'to go: ', line_to_go
            if heatmap_line.startswith(line_to_go):
                return True
        return False

    def format_heatmap(self):
        """
        We want to display the heatmap HTML in a Div
        Remove html, head, meta, and body tags
        """
        if self.has_error:
            return

        self.heatmap_html = '\n'.join(self.heatmap_html_lines)

        """
        formatted_lines = []
        for hline in self.heatmap_html_lines:
            print hline
            print '-' * 40
            if self.is_line_to_remove(hline.strip()):
                print 'SCRAP IT!!!'
                continue    # Skip this line!
            formatted_lines.append(hline)

        self.heatmap_html = '\n'.join(formatted_lines)
        """
