"""
Variables with expected file names from a successful run
of the dataset pipeline as well as prediction scripts

Dataset file diretory/
        - gentb_status_feedback.py
            (customized for each run, uses variables from this page)
        - /output/
            - result.json   (pipeline analyze result)
            - matrix.csv    (pipeline analyze result)
            - heatmap.html  (predict script result)
"""
RESULT_OUTPUT_DIRECTORY_NAME = 'output'    # local directory within the Dataset file directory
RESULT_JSON_FILE_NAME = 'result.json'
MATRIX_CSV_FILE_NAME = 'matrix.csv'
HEATMAP_HTML_FILE_NAME = 'matrix_heatmap.html'

RESULT_FILE_NAME_DICT = { 'RESULT_JSON_FILE_NAME' : RESULT_JSON_FILE_NAME,
            'MATRIX_CSV_FILE_NAME' : MATRIX_CSV_FILE_NAME,
            'HEATMAP_HTML_FILE_NAME' : HEATMAP_HTML_FILE_NAME,
            }

EXPECTED_FILE_NAME_LIST = (RESULT_JSON_FILE_NAME,
                            MATRIX_CSV_FILE_NAME,
                            HEATMAP_HTML_FILE_NAME)
EXPECTED_FILE_DESCRIPTIONS = ('results file',
                            'matrix file',
                            'heatmap file')
