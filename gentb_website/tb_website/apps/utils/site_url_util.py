from django.contrib.sites.models import Site

def get_site_url():

    current_site = Site.objects.get_current()
    return current_site.domain

