"""
Variables with expected file names from a successful run
of the dataset pipeline as well as prediction scripts

Dataset file diretory/
        - feedback.conf
            (customized for each run, uses variables from this page)
        - /output/
            - result.json   (predict script result)
            - matrix.csv    (pipeline analyze result)
            - heatmap.html  (predict script result)
"""
# local directory within the Dataset file directory
RESULT_OUTPUT_DIRECTORY_NAME = 'output'
RESULT_JSON_FILE_NAME = 'result.json'
MATRIX_JSON_FILE_NAME = 'matrix.json'

RESULT_FILE_NAME_DICT = {
    'RESULT_JSON_FILE_NAME' : RESULT_JSON_FILE_NAME,
    'MATRIX_JSON_FILE_NAME' : MATRIX_JSON_FILE_NAME,
}

EXPECTED_FILE_NAME_LIST = (
    RESULT_JSON_FILE_NAME,
    MATRIX_JSON_FILE_NAME,
)
EXPECTED_FILE_DESCRIPTIONS = [
    'results file',
    'matrix file',
]
