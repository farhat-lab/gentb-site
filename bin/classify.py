#
# Copyright (C) 2019  Dr. Maha Farhat
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# AUTHOR: Luca Freschi and David Hughes

import sys
import re


# Iterates through stdin looking for the string that matches kraken2's
# classification report. Output will look like this:
#   Loading database information... done.
#   1452549 sequences (680.14 Mbp) processed in 13.338s (6534.0 Kseq/m, 3059.49 Mbp/m).
#       1444500 sequences classified (99.45%) <---- matches this line
#       8049 sequences unclassified (0.55%)
#
# Then finds the percentage classified, and returns 0 only if the percentage
# classfied is greater than or equal to 90%.
def main():
    classified = re.compile(r'(?<!(un))classified \([0-9]+\.[0-9]*%\)')
    for line in sys.stdin:
        match = classified.search(line)
        if match:
            perc = re.compile(r'[0-9]+\.[0-9]*')
            percent = perc.search(match.group(0))
            tbperc = float(percent.group(0)) / 100

            if tbperc < 0.9:
                sys.exit(1)

            break


if __name__ == '__main__':
    main()
