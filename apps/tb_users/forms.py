from django.utils.translation import ugettext as _

from django import forms
from django.contrib.auth.models import User
from django.contrib.auth.forms import UserCreationForm

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


def get_signup_form_test_data():
    return dict(username='bmurray',
        password='hello123',
        retype_password='hello123',
        first_name='bill',
        last_name='murray',
        email='bill_murray@b.com',
        affiliation='iqss')
