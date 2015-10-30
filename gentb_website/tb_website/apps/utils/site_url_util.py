from django.conf import settings

from django.contrib.sites.models import Site

def get_site_url():

    current_site = Site.objects.get_current()

    if settings.IS_HTTPS_SITE:
        protocol = 'https://'
    else:
        protocol = 'http://'

    return protocol + current_site.domain
