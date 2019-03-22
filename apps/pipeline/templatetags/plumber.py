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


def replace_helper(text, fn, template, name):
    to_replace = r'(^| |=){}( |$)'.format(fn)
    return re.sub(to_replace, r'\1' + template.format(fn, name) + r'\2', text)


@register.filter("process_command")
def command(job):
    """processes the command by highlighting the input and output files"""
    text = job.debug_text
    input_fn = job.input_fn
    output_fn = job.output_fn

    # keep track of which files are inputs and which are outputs
    # set lookup has average case O(1)
    inputs = set()
    for fn in job.input_filenames():
        inputs.add(('input', fn))

    outputs = set()
    for fn in job.output_filenames():
        inputs.add(('output', fn))

    all_files = sorted(inputs | outputs, key=lambda tup: len(tup[1]))

    for io, fn in all_files:
        if io == 'input':
            name = fn.split('/files/')[-1].replace('XX:', '')
            if fn.startswith('XX:') or (job.pk and fn not in input_fn):
                text = text.replace(fn, INMIA.format(fn, name))
            elif '/bin/' in fn:
                name = name.split('/bin/')[-1]
                text = replace_helper(text, fn, INBIN, name)
            else:
                text = replace_helper(text, fn, INOK, name)
        else:
            name = fn.split('/')[-1]
            if fn not in output_fn and job.pk:
                text = replace_helper(text, fn, OUTMIA, name)
            else:
                text = replace_helper(text, fn, OUT, name)

    return mark_safe(text)


BINS = re.compile(r'(?P<io>\$){bin}(?P<suffix>[^\s;|>]*)')
INPUTS = re.compile(r'((?P<io>\$)(?P<prefix>[^{]*)' \
                    r'{(?P<literal>"?)(?P<name>[-\w]+)"?\}(?P<suffix>[^\s;|>]*))')
OUTPUTS = re.compile(r'((?P<io>\@)(?P<prefix>[^{]*)' \
                     r'{(?P<literal>"?)(?P<name>[-\w]+)"?\}(?P<suffix>[^\s;|>]*))')


@register.filter("process_template")
def template(program):
    """processes the template by highlighting the input and output files"""
    text = program.command_line
    text = BINS.sub(r'<span class="bin">\2</span>', text)
    text = INPUTS.sub(r'<span class="input">\1</span>', text)
    text = OUTPUTS.sub(r'<span class="output">\1</span>', text)
    return mark_safe(text)
