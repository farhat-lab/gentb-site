from django.conf.urls import patterns, include, url
from django.views.generic import TemplateView as Tv

urlpatterns = patterns('apps.basic_pages.views',
    url(r'^$', Tv.as_view(template_name='homepage.html'), name="view_homepage"),
    url(r'^about/$', Tv.as_view(template_name='about.html'), name="view_about_page"),
    url(r'^share/$', Tv.as_view(template_name='data_page.html'), name="view_data_page"),
    url(r'^terms-of-use/$', Tv.as_view(template_name='terms_of_use/terms-of-use-page.html'), name="view_terms_of_use"),
    url(r'^data-upload-information/$', Tv.as_view(template_name='data_upload_instructions.html'), name="view_data_upload_information"),

    url(r'^explore/$', 'view_explore_page', name="view_explore_page"),
)
