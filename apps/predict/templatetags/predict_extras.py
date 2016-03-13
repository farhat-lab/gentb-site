from django import template
register = template.Library()

@register.filter(name='readable_filesize')
def readable_filesize(value):
    """Convert bytes to a human readable string
    @param string or integer  file size in bytes
    @return string "human readable" filesize, e.g. 2.4 MB
            on fail: return 'n/a'
    """
    #print 'readable_filesize:', value
    num = value
    fail_str = 'n/a'
    if num.__class__.__name__ == 'str':
        if num.isdigit():
            num = int(num)
        else:
            return fail_str

    if num is None:
        #logger.error('Error: number for human readable filesize is "None"' % num)
        return 'n/a'

    for x in ['bytes','KB','MB','GB','TB']:
        try:
            if num < 1024.0:
                return "%3.1f %s" % (num, x)
            num /= 1024.0
        except:
            #logger.error('Error: could not convert %s to human readable filesize' % num)
            return fail_str

    return fail_str