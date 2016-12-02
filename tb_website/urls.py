"""
Global anchor for the website's urls
"""

from django.conf.urls import patterns, include, url
from django.conf.urls.static import static
from django.conf import settings
from django.views.generic import TemplateView as Tv

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    url(r'^$', Tv.as_view(template_name='home.html'), name="home"),
    url(r'^about/$', Tv.as_view(template_name='about.html'), name="about"),
    url(r'^share/$', Tv.as_view(template_name='share.html'), name="share"),
    url(r'^terms/$', Tv.as_view(template_name='terms.html'), name="terms"),
    url(r'^info/$', Tv.as_view(template_name='info.html'), name="info"),

    url(r'^explore/', include('apps.explore.urls', namespace='explore')),
    url(r'^predict/', include('apps.predict.urls', namespace='predict')),
    url(r'^genes/', include('apps.mutations.urls', namespace='genes')),
    url(r'^user/', include('apps.tb_users.urls', namespace='users')),
    url(r'^auth/', include('apps.tb_users.auth_urls')),
    url(r'^maps/', include('apps.maps.urls', namespace='maps')),
    url(r'^models/', include('django_spaghetti.urls')),

    url(r'^gentb-admin/', include(admin.site.urls)),
)

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

    import debug_toolbar
    urlpatterns += patterns('', url(r'^__debug__/', include(debug_toolbar.urls)))
