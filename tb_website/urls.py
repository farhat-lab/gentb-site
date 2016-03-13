from django.conf.urls import patterns, include, url
from django.conf.urls.static import static
from django.conf import settings
from django.views.generic import TemplateView

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    url(r'^', include('apps.basic_pages.urls')),
    url(r'^user/', include('apps.tb_users.urls')),
    url(r'^predict/', include('apps.predict.urls')),
    url(r'^maps/', include('apps.maps.urls')),

    url(r'^accounts/password_reset/$',
        'django.contrib.auth.views.password_reset',
        {'post_reset_redirect' : '/accounts/password_reset/mailed/'},
        name="password_reset"),

    (r'^accounts/password_reset/mailed/$',
        'django.contrib.auth.views.password_reset_done'),
        
    url(r'^accounts/password_reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
        'django.contrib.auth.views.password_reset_confirm',
        {'post_reset_redirect' : '/accounts/password_reset/complete/'},
                name="password_reset_confirm"),

    (r'^accounts/password_reset/complete/$',
        'django.contrib.auth.views.password_reset_complete'),

    url(r'^gentb-admin/', include(admin.site.urls)),
)

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

    import debug_toolbar
    urlpatterns += patterns('', url(r'^__debug__/', include(debug_toolbar.urls)))
