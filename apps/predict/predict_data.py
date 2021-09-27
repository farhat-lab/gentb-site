#
# Copyright (C) 2021 Maha Farhat
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
Because the json formatted output for prediction is a list and this list
has changed format several times (unfortunatly), this mask allows us to pick
out which data format we're dealing with.
"""

import json
from .utils import filter_none

def _float(val):
    if not val or str(val).lower() == 'none':
        return None
    if '.' in str(val):
        return float(val)
    raise ValueError("Not a float")

def _bool(val):
    if val in (True, False, 1, 0, '1', '0'):
        return bool(int(val))
    if str(val).lower() == 'true':
        return True
    if str(val).lower() == 'false':
        return False
    raise ValueError("Not a boolean")

# List[('name', type/test), ...]
PR_FORMATS = (
    [('drug_code', str), ('dr', float), ('fneg', float), ('fpos', float)],
    [('drug_code', str), ('pr', int), ('fneg', float), ('fpos', float), ('dr', float)],
)
if len([len(pr) for pr in PR_FORMATS]) != len(PR_FORMATS):
    raise ValueError("Prediction formats MUST be different lengths for detection")

class PredictParsingError(ValueError):
    pass

def decypher_predict_format(dat):
    metadata = {}
    rest, *data = dat
    for pr in PR_FORMATS:
        if len(rest) == len(pr):
            for (name, kind), datum in zip(pr, rest):
                metadata[name] = kind(datum)
            return metadata, data
    raise PredictParsingError("Can't decypher predict format!")


def prediction_from_file(matrix_fn):
    m_A, m_B, m_C, m_D = {}, {}, {}, {}
    with open(matrix_fn, 'r') as fhl:
        parts = json.loads(fhl.read())
        (pr, m_A, m_B) = parts[:3]

        # this is different because it assumes that there's only one strain.
        name = list(pr)[0][0]
        if not m_A:
            m_A[name] = {name: [[None] * len(pr)]}
        if not m_B:
            m_B[name] = {name: [[None] * len(pr)]}

        m_C[name] = [([None] * len(m_A[name][0]))] * len(m_A[name])
        m_D[name] = [([None] * len(m_A[name][0]))] * len(m_A[name])

        ex_1, ex_2 = parts[-2:] if len(parts) == 5 else ({}, {})
        for index, value in ex_1.items():
            if index < len(m_C[name]):
                m_C[name][int(index)] = filter_none(value)
        for index, value in ex_2.items():
            if index < len(m_D[name]):
                m_D[name][int(index)] = filter_none(value)

    m_A = list(zip(*m_A[name]))
    m_B = list(zip(*m_B[name]))
    m_C = list(zip(*m_C[name]))
    m_D = list(zip(*m_D[name]))
    # Rotate mutation matrix 90 degrees
    for x, (name, *rest) in enumerate(pr):
        yield (name, (list(rest), m_A[x], m_B[x], m_C[x], m_D[x]))
