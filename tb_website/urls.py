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
    url(r'^data/$', Tv.as_view(template_name='data.html'), name="data"),
    url(r'^data/info/$', Tv.as_view(template_name='info.html'), name="info"),
    url(r'^terms/$', Tv.as_view(template_name='terms.html'), name="terms"),

    url(r'^gentb-admin/', include(admin.site.urls)),
    url(r'^models/', include('django_spaghetti.urls', namespace='spaghetti')),
    url(r'^user/', include('apps.tb_users.urls', namespace='users')),
    url(r'^auth/', include('apps.tb_users.auth_urls')),

    #url(r'.+', Tv.as_view(template_name='offline.html'), name="offline"),

    url(r'^explore/', include('apps.explore.urls', namespace='explore')),
    url(r'^predict/', include('apps.predict.urls', namespace='predict')),
    url(r'^pipeline/', include('apps.pipeline.urls', namespace='pipeline')),
    url(r'^uploads/', include('apps.uploads.urls', namespace='uploads')),
    url(r'^genes/', include('apps.mutations.urls', namespace='genes')),
    url(r'^maps/', include('apps.maps.urls', namespace='maps')),

)

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

    import debug_toolbar
    urlpatterns += patterns('', url(r'^__debug__/', include(debug_toolbar.urls)))


from .views import Error

for e in ('403','404','500'):
    locals()['handler'+e] = Error.as_error(e)
    urlpatterns += patterns('', url('^error/%s/$' % e, Error.as_error(e)))

