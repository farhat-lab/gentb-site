"""
Variables with expected file names from a successful run
of the dataset pipeline as well as prediction scripts

Dataset file diretory/
        - /output/
            - matrix.csv    (pipeline analyze result)
            - matrix.json   (heatmap related info)
"""
# local directory within the Dataset file directory
RESULT_OUTPUT_DIRECTORY_NAME = 'output'
MATRIX_JSON_FILE_NAME = 'matrix.json'

RESULT_FILE_NAME_DICT = {
    'MATRIX_JSON_FILE_NAME' : MATRIX_JSON_FILE_NAME,
}

EXPECTED_FILE_NAME_LIST = (
    MATRIX_JSON_FILE_NAME,
)
EXPECTED_FILE_DESCRIPTIONS = [
    'matrix file',
]
