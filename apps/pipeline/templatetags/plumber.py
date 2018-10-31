"""
Provide processing for pipeline views
"""
import re

from django.utils.safestring import mark_safe
from django.template import Library

register = Library() # pylint: disable=invalid-name

INMIA = '<span class="input-missing" title="Missing input: {}" data-toggle="tooltip" data-placement="bottom">{}</span>'
INOK = '<span class="input" title="{}" data-toggle="tooltip" data-placement="bottom">{}</span>'
INBIN = '<span class="bin" title="{}" data-toggle="tooltip" data-placement="bottom">{}</span>'
OUT = '<span class="output" data-toggle="tooltip" data-placement="bottom" title="{}">{}</span>'
OUTMIA = '<span class="input-missing" data-toggle="tooltip" data-placement="bottom" title="Missing output: {}">{}</span>'

@register.filter("process_command")
def command(job):
    """processes the command by highlighting the input and output files"""
    text = job.debug_text
    input_fn = job.input_fn
    output_fn = job.output_fn

    for ifn in set(job.input_filenames()):
        name = ifn.split('/files/')[-1].replace('XX:', '')
        if ifn.startswith('XX:') or (job.pk and ifn not in input_fn):
            text = text.replace(ifn, INMIA.format(ifn, name))
        elif '/bin/' in ifn:
            name = name.split('/bin/')[-1]
            text = text.replace(ifn, INBIN.format(ifn, name))
        else:
            text = text.replace(ifn, INOK.format(ifn, name))

    for ofn in set(job.output_filenames()):
        name = ofn.split('/')[-1]
        if ofn not in output_fn and job.pk:
            text = text.replace(ofn, OUTMIA.format(ofn, name))
        else:
            text = text.replace(ofn, OUT.format(ofn, name))

    return mark_safe(text)

BINS = re.compile(r'(?P<io>\$){bin}(?P<suffix>[^\s;|>]*)')
INPUTS = re.compile(r'((?P<io>\$)(?P<prefix>[^{]*)' \
                    r'{(?P<name>[-\w]+)\}(?P<suffix>[^\s;|>]*))')
OUTPUTS = re.compile(r'((?P<io>\@)(?P<prefix>[^{]*)' \
                    r'{(?P<name>[-\w]+)\}(?P<suffix>[^\s;|>]*))')

@register.filter("process_template")
def template(program):
    """processes the template by highlighting the input and output files"""
    text = program.command_line
    text = BINS.sub(r'<span class="bin">\2</span>', text)
    text = INPUTS.sub(r'<span class="input">\1</span>', text)
    text = OUTPUTS.sub(r'<span class="output">\1</span>', text)
    return mark_safe(text)
