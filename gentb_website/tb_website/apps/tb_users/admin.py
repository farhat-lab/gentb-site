from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from django.contrib.auth.models import User

from apps.tb_users.models import TBUser, TBAdminContact

# Define an inline admin descriptor for Employee model
# which acts a bit like a singleton
class TBUserInline(admin.StackedInline):
    model = TBUser
    can_delete = False
    verbose_name_plural = 'TB user'

# Define a new User admin
class UserAdmin(UserAdmin):
    inlines = (TBUserInline, )


class TBUserAdmin(admin.ModelAdmin):
    save_on_top = True
    #search_fields = ('title', 'affiliation', 'last_name', 'first_name')
    #list_filter = [ 'affiliation', 'has_prediction']
    readonly_fields = [ 'md5', ]
admin.site.register(TBUser, TBUserAdmin)

class TBAdminContactAdmin(admin.ModelAdmin):
    save_on_top = True
    list_display = ('email', 'name', 'is_active')
    #search_fields = ('title', 'affiliation', 'last_name', 'first_name')
    #list_filter = [ 'affiliation', 'has_prediction']
    readonly_fields = [ 'created', 'modified' ]
admin.site.register(TBAdminContact, TBAdminContactAdmin)


# Re-register UserAdmin to include extra TBUser fields
#
admin.site.unregister(User)
admin.site.register(User, UserAdmin)
