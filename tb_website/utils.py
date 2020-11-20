
from math import log

UNITS = ['bytes', 'kB', 'MB', 'GB', 'TB', 'PB']
DECIMS = [0, 0, 1, 2, 2, 2]

def sizeof(num):
    if num > 1:
        exponent = min(int(log(num, 1024)), len(UNITS) - 1)
        quotient = float(num) / 1024 ** exponent
        unit = UNITS[exponent]
        num_decimals = DECIMS[exponent]
        format_string = '{:.%sf} {}' % (num_decimals)
        return format_string.format(quotient, unit)
    if num == 0:
        return '0 bytes'
    if num == 1:
        return '1 byte'

def to(method):
    """Turn generators into objects, method can be a type, obj or function"""
    def _outer(f):
        def _inner(*args, **kw):
            return method(f(*args, **kw))
        return _inner
    return _outer


