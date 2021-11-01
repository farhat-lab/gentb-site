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


# trims a VCF file by excluding any rows where the ALT is '.' or the REF has
# only one base
def main():
    if len(sys.argv) < 2:
        sys.exit(1)

    REF_IDX = -1
    ALT_IDX = -1

    f = open(sys.argv[1], 'r')

    for line in f:
        if line.startswith('#'):
            sys.stdout.write(line)

            if line.startswith('##'):
                continue

            # extract which column REF and ALT are in
            headers = line.rstrip('\n').split('\t')
            REF_IDX = headers.index('REF')
            ALT_IDX = headers.index('ALT')

            continue

        data = line.rstrip('\n').split('\t')

        # if REF_IDX and ALT_IDX aren't set then this will fail
        if len(data[REF_IDX]) != 1 or data[ALT_IDX] != '.':
            sys.stdout.write(line)

    f.close()


if __name__ == '__main__':
    main()
