from django.conf import settings

from django.contrib.sites.models import Site

def get_site_url(for_internal_callback=False):

    current_site = Site.objects.get_current()

    # For HMS prod, an internal callback is used
    #
    if for_internal_callback:
        if settings.INTERNAL_CALLBACK_SITE_URL:
            return settings.INTERNAL_CALLBACK_SITE_URL
            #return 'http://rc-app-shared01.orchestra:9001'

    if settings.IS_HTTPS_SITE:
        protocol = 'https://'
    else:
        protocol = 'http://'

    return protocol + current_site.domain
