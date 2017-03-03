#!/usr/bin/env python
#
# Copyright (C) 2017 Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Quality control script, will exit with a non-zero status if any qc file
or stdin file results in a failure.

Prints error information to stderr and echos any stdin qc information
back to stdout for pipeline processing.
"""

import os
import sys
import atexit


def test_quality(name, fhl, threshold=15):
    """
    Parse the quality informaion from samtools output
    """
    sys.stderr.write(name)
    sys.stderr.write("\t")

    # Get each of the quality values from the qc format
    qc = list([int(line.rstrip().split('\t')[2]) for line in fhl if '\t' in line])

    if not qc:
        sys.stderr.write('QC file empty!\n')
        return 1

    # Get the percentage of quality values below the threshold threshold
    pc = (len([1 for depth in qc if depth < threshold]) * 100) / len(qc)

    if pc < 2:
        sys.stderr.write('PASS\n')
        return 0
    else:
        sys.stderr.write('FAIL : %d%% of the DR sites have coverage <%dx\n' % (pc, threshold))
        return 1


def test_qc_file(filename, threshold=15):
    """
    Parse a qc file saved from the samtools process.
    """
    name = os.path.basename(filename).rsplit('.', 1)[0]
    with open(filename, 'r') as fhl:
        return test_quality(name, fhl, threshold=threshold)
        

if __name__ == '__main__':
    threshold = 15
    failed = 0
    for arg in sys.argv[1:]:
        if arg.startswith('-Q='):
            threshold = int(arg.split('=')[-1])
            continue
        elif os.path.isdir(arg):
            for filename in os.listdir(arg):
                if filename.endswith('.qc'):
                    failed += test_qc_file(filename, threshold)
        elif os.path.isfile(arg):
            failed += test_qc_file(arg, threshold)
        else:
            sys.stderr.write("File not found: '%s'\n" % arg)
            sys.exit(1)

    # Detect piped in data
    if not sys.stdin.isatty():
        echo = list(sys.stdin)
        sys.stdout.write("".join(echo))
        failed += test_quality("stdin", echo, threshold)

    if failed:
        sys.stderr.write("FAILURES: %d\n" % failed)
        sys.exit(5)


