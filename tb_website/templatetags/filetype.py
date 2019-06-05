"""
File type tools
"""

from collections import OrderedDict
from django.template import Library

register = Library() # pylint: disable=invalid-name

TYPES = OrderedDict([
    ('vcf', ('vcf', 'evcf')),
    ('sam', ('sam',)),
    ('bam', ('bam', 'bai')),
    ('fastq', ('fastq', 'fasta')),

    ('spreadsheet', ('ods', 'xls', 'xlsx', 'csv', 'tsv')),
    ('archive', ('zip', 'tar', 'gz', 'rar', 'xz', 'bz2')),
    ('image', ('png', 'gif', 'jpg', 'jpeg', 'tiff')),
    ('audio', ('wav', 'mp3', 'aiff')),
    ('code', ('py', 'pl', 'r', 'json')),
    ('html', ('html', 'htm')),
    ('text', ('txt', 'var')),
    ('pdf', ('pdf',)),
    ('xml', ('xml',)),
])

@register.filter("fileicon")
def get_fileicon(filename):
    """Returns the static location of of the file icon for this type"""
    if filename.endswith('.gz'):
        filename = filename[:-3]
    ext = filename.rsplit('.', 1)[-1]
    for (icon, lst) in TYPES.items():
        if ext in lst:
            return 'images/files/{}.svg'.format(icon)
    return 'images/files/file.svg'
