from apps.predict.models import PredictDataset
from django.views.decorators.cache import cache_page
from apps.utils.site_url_util import get_site_url

USERNAME_KEY = 'username'
IS_LOGGED_IN_KEY = 'IS_LOGGED_IN'
IS_ACTIVE_STAFF = 'IS_ACTIVE_STAFF'
GENTB_CONTACT_EMAIL = 'TBpredict@gmail.com'
GENTB_CONTACT_EMAIL_HTML = '<a href="mailto:{0}">{0}</a>'.format(GENTB_CONTACT_EMAIL)

def get_common_dict(request, title, **kwargs):


    d = dict(title=title,
             IS_LOGGED_IN=False,
             SITE_URL=get_site_url(),
             GENTB_CONTACT_EMAIL_HTML=GENTB_CONTACT_EMAIL_HTML,
             GENTB_CONTACT_EMAIL=GENTB_CONTACT_EMAIL)

    # Add other arguments
    for k, v in kwargs.items():
        d[k] = v

    if request.user.is_authenticated():
        d[USERNAME_KEY] = request.user.username
        d[IS_LOGGED_IN_KEY] = True
        d['my_dataset_count'] = get_dataset_count(request)
        d[IS_ACTIVE_STAFF] = request.user.is_active and request.user.is_staff

    return d



def get_dataset_count(request):

    if not request.user.is_authenticated():
        return 0

    if not hasattr(request.user, 'tbuser'):
        return 0

    return PredictDataset.objects.filter(user=request.user.tbuser).count()
