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
# AUTHOR: David Hughes

import sys


# check what proportion of depths are >= 10
# want at least 95% >= 10
def main():
    num_lines = 0
    valid = 0

    for line in sys.stdin:
        num_lines += 1

        # samtools depth -a format has depth in third column
        depth = int(line.split('\t')[2])
        if depth >= 10:
            valid += 1

    # nonzero exit for GenTB Pipeline
    if (valid / num_lines) < 0.95:
        sys.exit(1)


if __name__ == '__main__':
    main()
