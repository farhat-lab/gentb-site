from __future__ import print_function
from datetime import datetime, timedelta
from django.utils import timezone

from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext

from apps.tb_users.models import TBUser
from apps.tb_users.forms import SignUpForm, get_signup_form_test_data
from apps.utils.view_util import get_common_dict



def view_signup_page(request):

    d = get_common_dict(request, 'Sign Up', signup_page=True)

    # If the person is already logged in, go to the homepage
    #
    if request.user.is_authenticated():
        return HttpResponseRedirect(reverse('view_homepage', kwargs={}))

    if request.POST:
        f = SignUpForm(request.POST)
        if f.is_valid():
            tb_user = f.create_tb_user()
            #return HttpResponse('TB user created: %s' % tb_user)
            success_url = reverse('view_signup_success',
                                 kwargs=dict(tb_user_md5=tb_user.md5)
                                )
            return HttpResponseRedirect(success_url)
        else:
            d['ERROR_FOUND'] = True
    else:
        #test_data = get_signup_form_test_data()
        #f = SignUpForm(initial=test_data)
        f = SignUpForm()

    d['signup_form'] = f

    return render_to_response('tb_users/signup_page.html',
                             d,
                             context_instance=RequestContext(request))


def view_signup_success(request, tb_user_md5):


    d = get_common_dict(request, 'Sign Up', signup_page=True)

    try:
        tb_user = TBUser.objects.get(md5=tb_user_md5)
    except TBUser.DoesNotExist:
        raise Http404('User not found')

    # This page should only be seen once, expire it in 10 minutes
    #
    right_now = timezone.now()  #datetime.now()
    elapsed = right_now - tb_user.user.date_joined
    if elapsed > timedelta(minutes=15):
        raise Http404('Page not found')

    d['tb_user'] = tb_user.user
    d['tb_user_affiliation'] = tb_user.affiliation

    return render_to_response('tb_users/signup_success_page.html',
                             d,
                             context_instance=RequestContext(request))
