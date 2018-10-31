"""
Provide processing for pipeline views
"""
import re

from django.utils.safestring import mark_safe
from django.template import Library

register = Library() # pylint: disable=invalid-name

INMIA = '<span class="input-missing" title="Missing file: {}" data-toggle="tooltip" data-placement="bottom">{}</span>'
INOK = '<span class="input" title="{}" data-toggle="tooltip" data-placement="bottom">{}</span>'
INBIN = '<span class="bin" title="{}" data-toggle="tooltip" data-placement="bottom">{}</span>'
OUT = '<span class="output" data-toggle="tooltip" data-placement="bottom" title="{}">{}</span>'

@register.filter("process_command")
def command(run):
    """processes the command by highlighting the input and output files"""
    text = run.debug_text
    for ifn in set(run.input_filenames()):
        name = ifn.split('/files/')[-1].replace('XX:', '')
        if ifn.startswith('XX:'):
            text = text.replace(ifn, INMIA.format(ifn, name))
        elif '/bin/' in ifn:
            name = name.split('/bin/')[-1]
            text = text.replace(ifn, INBIN.format(ifn, name))
        else:
            text = text.replace(ifn, INOK.format(ifn, name))

    for ofn in set(run.output_filenames()):
        name = ofn.split('/')[-1]
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
