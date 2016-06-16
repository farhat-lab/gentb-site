
from django.views.generic import CreateView, TemplateView
from django.core.urlresolvers import reverse_lazy
from apps.tb_users.forms import SignUpForm

class UserSignUp(CreateView):
    template_name = 'tb_users/signup_page.html'
    success_url = reverse_lazy('users:signup_success')
    form_class = SignUpForm


class SignUpSuccess(TemplateView):
    template_name = 'tb_users/signup_success_page.html'

