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

    #url(r'^share/', include('apps.shared_data.urls')),

    url(r'^predict/', include('apps.predict.urls')),

    url(r'^maps/', include('apps.maps.urls')),

    # Examples:
    # url(r'^$', 'tb_website.views.home', name='home'),
    # url(r'^tb_website/', include('tb_website.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^tb-admin/', include(admin.site.urls)),
)

# Uncomment the next line to serve media files in dev.
urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

if settings.DEBUG:
    import debug_toolbar
    urlpatterns += patterns('',
                            url(r'^__debug__/', include(debug_toolbar.urls)),
                            )
