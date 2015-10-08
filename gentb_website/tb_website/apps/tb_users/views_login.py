from __future__ import print_function

from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.template import RequestContext

from django.contrib.auth import login, logout

from apps.tb_users.forms import TBUserAuthenticationForm

from apps.utils.view_util import get_common_dict, IS_LOGGED_IN_KEY
from apps.utils.msg_util import *

from apps.basic_pages.views import view_homepage

def view_homepage_just_logged_in(request):

    return view_homepage(request, just_logged_in=True)


def view_login_after_logout(request):

    return view_login_page(request, just_logged_out=True)

def view_login_page(request, just_logged_out=False):

    d = get_common_dict(request, 'Log In', login_page=True)
    d['JUST_LOGGED_OUT'] = just_logged_out

    # If the person is already logged in, go to the homepage
    #
    if request.user.is_authenticated():
        d['ALREADY_LOGGED_IN'] = True
        return render_to_response('tb_users/login_page.html',
                             d,
                             context_instance=RequestContext(request))

    if request.POST:
        f = TBUserAuthenticationForm(data=request.POST)
        if f.is_valid():
            user_obj = f.get_user()

            # Log the person in
            login(request, f.get_user())

            # Go to log in success
            success_url = reverse('view_homepage_just_logged_in', kwargs={})
            return HttpResponseRedirect(success_url)

        else:
            d['ERROR_FOUND'] = True
    else:
        f = TBUserAuthenticationForm()

    d['login_form'] = f

    return render_to_response('tb_users/login_page.html',
                             d,
                             context_instance=RequestContext(request))



def view_logout_page(request):

    d = get_common_dict(request, 'Log Out', logout_page=True)

    # Not logged in, redirect to the log in page
    if d[IS_LOGGED_IN_KEY] is False:
        d['NOT_LOGGED_IN'] = True
        login_url = reverse('view_login_page',
                                 kwargs={})

        return HttpResponseRedirect(login_url)

    # Log the person out
    logout(request)

    return view_login_after_logout(request)



