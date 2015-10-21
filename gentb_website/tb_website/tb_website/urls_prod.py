from django.conf.urls import patterns, include, url
from django.conf.urls.static import static
from django.conf import settings
from django.views.generic import TemplateView

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
#url(r'^$', TemplateView.as_view(template_name='base.html')),

  url(r'^', include('apps.basic_pages.urls')),

  url(r'^user/', include('apps.tb_users.urls')),

  url(r'^predict/', include('apps.predict.urls')),

  url(r'^maps/', include('apps.maps.urls')),

  # Examples:
  # url(r'^$', 'tb_website.views.home', name='home'),
  # url(r'^tb_website/', include('tb_website.foo.urls')),

  # Uncomment the admin/doc line below to enable admin documentation:
  # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

  # Uncomment the next line to enable the admin:
  url(r'^gentb-admin/', include(admin.site.urls)),
)
