"""
Accepted File patterns for GenTB.  Used with dropbox links.
"""
import re

# Used by apps/predict/models.py
#
FILE_TYPE_VCF = 'vcf'
FILE_TYPE_FASTQ = 'fastq'
FILE_TYPES = [(FILE_TYPE_VCF, 'VCF'), (FILE_TYPE_FASTQ, 'FastQ')]

FASTQ_PAIR_ENDED = 'pair-end'
FASTQ_SINGLE_ENDED = 'single-end'
FASTQ_FILE_TYPES = [(FASTQ_PAIR_ENDED, 'Pair-end'),\
    (FASTQ_SINGLE_ENDED, 'Single-end')]

FASTQ_PAIR_END_EXTENSION_R = '_R'
FASTQ_PAIR_END_EXTENSION_DOT = '.'


# e.g. teststrain_R1.fastq & teststrain_R2.fastq
FASTQ_PAIR_END_EXTENSION_R_PATTERNS = [r'_R\d{1,9}\.fastq$', r'_R\d{1,9}\.fastq\.']

# e.g. teststrain.1.fastq & teststarin.2.fastq
FASTQ_PAIR_END_EXTENSION_DOT_PATTERNS = [r'\.\d{1,9}\.fastq$', r'\.\d{1,9}\.fastq\.']

FASTQ_PAIR_END_EXTENSION_TYPES = [(FASTQ_PAIR_END_EXTENSION_R, FASTQ_PAIR_END_EXTENSION_R),\
                    (FASTQ_PAIR_END_EXTENSION_DOT, FASTQ_PAIR_END_EXTENSION_DOT)]

GENTB_FASTQ_FILE_PATTERNS = [r'\.fastq$', r'\.fastq\.']
GENTB_VCF_FILE_PATTERNS = [r'\.vcf$', r'\.vcf\.']
GENTB_FILE_PATTERNS = GENTB_FASTQ_FILE_PATTERNS + GENTB_VCF_FILE_PATTERNS

class FilePatternHelper(object):


    @staticmethod
    def is_vcf_file(file_type):
        if file_type == FILE_TYPE_VCF:
            return True
        return False

    @staticmethod
    def is_fastq_file(file_type):
        if file_type == FILE_TYPE_FASTQ:
            return True
        return False

    @staticmethod
    def is_fastq_single_ended(fastq_type):
        if fastq_type == FASTQ_SINGLE_ENDED:
            return True
        return False

    @staticmethod
    def is_fastq_pair_ended(fastq_type):
        if fastq_type == FASTQ_PAIR_ENDED:
            return True
        return False


    @staticmethod
    def get_fastq_extension_type(list_of_filenames):
        """
        For pair-ended FastQ files, figure out
        if the extension type is "_R" or "."
        """
        if not list_of_filenames or len(list_of_filenames) == 0:
            return None

        # First make sure there are some FastQ files
        #
        found_fastq = False
        for fname in list_of_filenames:
            for pat in GENTB_FASTQ_FILE_PATTERNS:
                if re.search(pat, fname, re.IGNORECASE):
                    found_fastq = True

        if found_fastq is False:
            return None


        # Now look for the _R pattern as in:
        #   - teststrain_R1.fastq
        #   - teststrain_R2.fastq
        #
        for fname in list_of_filenames:
            for pat in FASTQ_PAIR_END_EXTENSION_R_PATTERNS:
                if re.search(pat, fname, re.IGNORECASE):
                    # Found it, return R pattern
                    return FASTQ_PAIR_END_EXTENSION_R

        # Default to "." pattern
        #
        #  Example: teststrain.1.fastq, teststarin.2.fastq
        for fname in list_of_filenames:
            for pat in FASTQ_PAIR_END_EXTENSION_DOT_PATTERNS:
                if re.search(pat, fname, re.IGNORECASE):
                    # Found it, return R pattern
                    return FASTQ_PAIR_END_EXTENSION_DOT

        return None


    @staticmethod
    def get_file_patterns_for_dropbox(file_type):
        """
        Based on the filetype, return regex search patterns used to
        check the dropbox link metadata
        """
        if file_type is None:
            return None

        if file_type == FILE_TYPE_VCF:
            return GENTB_VCF_FILE_PATTERNS
        elif file_type == FILE_TYPE_FASTQ:
            return GENTB_FASTQ_FILE_PATTERNS
        else:
            return None


    @staticmethod
    def get_file_patterns_err_msg(file_patterns):

        if file_patterns is None:
            return 'Please make sure you have at least one ".fastq" or ".vcf" file.'
        elif file_patterns == GENTB_FASTQ_FILE_PATTERNS:
            return 'Please make sure you have at least one ".fastq" file.'
        elif file_patterns == GENTB_VCF_FILE_PATTERNS:
            return 'Please make sure you have at least one ".vcf" file.'
        else:
            return 'Please make sure you have at least one ".fastq" or ".vcf" file.'
