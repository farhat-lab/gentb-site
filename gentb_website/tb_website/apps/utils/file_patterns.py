"""
Accepted File patterns for GenTB.  Used with dropbox links.
"""

# Used by apps/predict/models.py
#
FILE_TYPE_VCF = 'vcf'
FILE_TYPE_FASTQ = 'fastq'
FILE_TYPES = [(FILE_TYPE_VCF, 'VCF'), (FILE_TYPE_FASTQ, 'FastQ')]

FASTQ_PAIR_ENDED = 'pair-end'
FASTQ_SINGLE_ENDED = 'single-end'
FASTQ_FILE_TYPES = [(FASTQ_PAIR_ENDED, 'Pair-end'),\
    (FASTQ_SINGLE_ENDED, 'Single-end')]


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
