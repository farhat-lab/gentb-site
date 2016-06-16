
from os.path import join

from django.conf import settings
from django.contrib.sites.models import Site

VCF_ANALYSIS_SCRIPT = 'analyseVCF.pl'
FASTQ_ANALYSIS_SCRIPT = 'analyseNGS.pl'
MANUAL_ANALYSIS_SCRIPT = 'analyseMan.pl'
SCRIPT_DIR = join(settings.SITE_ROOT, 'apps', 'predict', 'pipeline')

FILE_TYPE_VCF = 'vcf'
FILE_TYPE_FASTQ = 'fastq'
FILE_TYPE_MANUAL = 'manual'
FILE_TYPES = [ 
  (FILE_TYPE_VCF, 'Variant Call Format (VCF)'),
  (FILE_TYPE_FASTQ, 'FastQ Nucleotide Sequence'),
  (FILE_TYPE_MANUAL, 'Mutations Manual Entry'),
]

FASTQ_PAIR_ENDED = 'pair-end'
FASTQ_SINGLE_ENDED = 'single-end'
FASTQ_FILE_TYPES = [ 
  (FASTQ_PAIR_ENDED, 'Pair-end'),
  (FASTQ_SINGLE_ENDED, 'Single-end'),
]
FASTQ_PAIR_END = { 
  FASTQ_PAIR_ENDED: '_R',
  FASTQ_SINGLE_ENDED: '.',
}

def get_site_url(internal=False):
    """
    Returns the right server address for this website.
    """
    protocol = 'http%s://' % ('', 's')[settings.IS_HTTPS_SITE]

    if internal and settings.INTERNAL_CALLBACK_SITE_URL:
        if '://' not in settings.INTERNAL_CALLBACK_SITE_URL:
            return protocol + settings.INTERNAL_CALLBACK_SITE_URL
        return settings.INTERNAL_CALLBACK_SITE_URL
    return protocol + Site.objects.get_current().domain

