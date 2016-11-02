from django.utils.translation import ugettext as _
import re

from django import forms
from django.contrib.auth.models import User
from django.contrib.auth.forms import AuthenticationForm, UserCreationForm

from apps.tb_users.models import TBUser

MINIMUM_PASSWORD_LENGTH = 7


class SignUpForm(UserCreationForm):

    error_messages = {
        'password_mismatch': _("The two password fields didn't match."),
        'password_complexity_length': _("The password must be at least 7 characters long (and contain at least 1 letter and 1 digit.)"),
        'password_complexity_content': _("The password must contain at least 1 letter and 1 digit."),
    }

    affiliation = forms.CharField(label='Affiliation', max_length=255)

    def save(self, **kw):
        # save User, made active
        user = super(SignUpForm, self).save(commit=False)
        user.is_active = True
        user.save()

        tb_user = TBUser(user=user, affiliation=self.cleaned_data['affiliation'])
        tb_user.save()

        return tb_user

    class Meta:
        model = User
        fields = ['username', 'password1', 'password2', 'first_name', 'last_name', 'email', 'affiliation']


class LoginForm(AuthenticationForm):

    def confirm_login_allowed(self, user):
        """
        Controls whether the given User may log in. This is a policy setting,
        independent of end-user authentication. This default behavior is to
        allow login by active users, and reject login by inactive users.
        If the given user cannot log in, this method should raise a
        ``forms.ValidationError``.
        If the given user may log in, this method should return None.
        """
        if not user.is_active:
            raise forms.ValidationError(
                self.error_messages['inactive'],
                code='inactive',
            )

        if not hasattr(user, 'tbuser'):
            raise forms.ValidationError(
                self.error_messages['not_tbuser'],
                code='not_tbuser',
            )

    def __init__(self, *args, **kwargs):
        super(LoginForm, self).__init__(*args, **kwargs)

        # Add additional error message
        self.error_messages['not_tbuser'] = 'Sorry! You are not a TBUser. '\
                                            'Please contact the administrator'


def get_signup_form_test_data():
    return dict(username='bmurray',
        password='hello123',
        retype_password='hello123',
        first_name='bill',
        last_name='murray',
        email='bill_murray@b.com',
        affiliation='iqss')
