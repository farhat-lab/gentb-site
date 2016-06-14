from django.conf import settings

from django.contrib.sites.models import Site

def get_site_url(internal=False):
    """
    Returns the right server address for this website.
    """
    if internal and settings.INTERNAL_CALLBACK_SITE_URL:
        return settings.INTERNAL_CALLBACK_SITE_URL

    protocol = 'http%s://' % ('', 's')[settings.IS_HTTPS_SITE]
    return protocol + Site.objects.get_current().domain
