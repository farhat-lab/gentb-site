#!/usr/bin/env python3
# pylint: disable=wrong-import-position
"""
Check all mutation names and format them as needed.
"""

import sys

sys.path.insert(0, '.')
sys.path.insert(0, '..')

try:
    import manage # pylint: disable=unused-import
except ImportError as err:
    sys.stderr.write("Could not run script! Is manage.py not in the current"\
        "working directory, or is the environment not configured?:\n"\
        "{:s}\n".format(err))
    sys.exit(1)

from argparse import ArgumentParser
from apps.mutations.models import Mutation

class Command():
    """Check all mutation names and format the data as needed."""
    def __init__(self):
        self.arg_parser = ArgumentParser(description=self.__doc__)
        self.add_arguments(self.arg_parser)
        self.args = self.arg_parser.parse_args(sys.argv[1:])

    def debug(self, msg):
        """Conditional debug message"""
        if self.args.debug:
            sys.stderr.write(msg + "\n")

    def add_arguments(self, pars):
        """Add extra command arguments"""
        pars.add_argument("--debug", action='store_true',\
            help='Print out more information about errors and what is happening.')
        return self

    def handle(self):
        """Handle the command being called"""
        unmatched = []
        for mutation in Mutation.objects.all().order_by('nucleotide_position', 'name'):
            self.debug(f" * Looking at {mutation.name}")
            try:
                changes = mutation.name_to_data()
                for m_field, original, data in changes:
                    self.debug(f"  > {m_field} changed '{original}' to '{data}'")
                if changes:
                    self.debug(f" ++ Saving mutation changes")
                    mutation.save()
            except ValueError:
                unmatched.append(mutation.name)
                continue
            except KeyError as err:
                self.debug(err)

        for mut in unmatched:
            sys.stderr.write(f" ! {mut}\n")


if __name__ == '__main__':
    Command().handle()
