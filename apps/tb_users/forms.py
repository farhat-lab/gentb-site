from django.utils.translation import ugettext as _
import re

from django import forms
from django.contrib.auth.models import User
from django.contrib.auth.forms import AuthenticationForm

from apps.tb_users.models import TBUser

MINIMUM_PASSWORD_LENGTH = 7

class SignUpForm(forms.ModelForm):

    error_messages = {
        'password_mismatch': _("The two password fields didn't match."),
        'password_complexity_length': _("The password must be at least 7 characters long (and contain at least 1 letter and 1 digit.)"),
        'password_complexity_content': _("The password must contain at least 1 letter and 1 digit."),
    }

    affiliation = forms.CharField(label='Affiliation', max_length=255)
    retype_password = forms.CharField(label='Retype Password', widget=forms.PasswordInput())

    def clean_email(self):
        email = self.cleaned_data.get('email')
        username = self.cleaned_data.get('username')
        if email and User.objects.filter(email=email).exclude(username=username).count():
            raise forms.ValidationError(u'Email addresses must be unique.')
        return email

    def clean_password(self):
        password = self.cleaned_data.get('password')

        if len(password) < MINIMUM_PASSWORD_LENGTH:
            raise forms.ValidationError(
                self.error_messages['password_complexity_length'],
                code='password_complexity_length',
        )

        #check if contains a digit
        if not re.search(r'\d', password):
            raise forms.ValidationError(
                self.error_messages['password_complexity_content'],
                code='password_complexity_content',
            )

        # check for at least 1 letter
        if not (re.search(r'[A-Z]', password) or re.search(r'[a-z]', password)):
            raise forms.ValidationError(
                self.error_messages['password_complexity_content'],
                code='password_complexity_content',
            )

        return password

    def clean_retype_password(self):
        password = self.cleaned_data.get('password')
        retype_password = self.cleaned_data.get('retype_password')
        if password and retype_password:
            if password != retype_password:
                raise forms.ValidationError(
                    self.error_messages['password_mismatch'],
                    code='password_mismatch',
                )
        return retype_password

    def create_tb_user(self):        
        assert hasattr(self, 'cleaned_data'), "Do not call this method on an invalid form."

        # save User, made inactive
        #
        user = self.save(commit=False)  # get object
        #user.is_active = False  # set to inactive, use for activation link scenario
        user.save() # save

        # save TBUser
        #
        tb_user = TBUser(user=user, affiliation=self.cleaned_data['affiliation'])
        tb_user.save()

        return tb_user

    class Meta:
        model = User
        fields = ['username', 'password', 'retype_password', 'first_name', 'last_name', 'email', 'affiliation']

    def __init__(self, *args, **kwargs):
        super(SignUpForm, self).__init__(*args, **kwargs)

        for key in self.fields:
            self.fields[key].required = True

        self.fields['password'].widget = forms.PasswordInput()
        #password = forms.CharField(widget=forms.PasswordInput())


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
