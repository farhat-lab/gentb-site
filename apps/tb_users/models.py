from datetime import datetime
from hashlib import md5

from django.db import models
from django.db.models import CASCADE
from django.contrib.auth.models import User

from model_utils.models import TimeStampedModel

class TBUser(models.Model):
    user = models.OneToOneField(User, related_name='tbuser', on_delete=CASCADE)
    affiliation = models.CharField(max_length=255)

    md5 = models.CharField(max_length=40, blank=True, db_index=True, help_text='auto-filled on save')

    def __str__(self):
        if self.user.last_name and self.user.first_name:
            return '%s, %s' % (self.user.last_name, self.user.first_name)
        return self.user.username

    def save(self, *args, **kwargs):
        if not self.id:
            super(TBUser, self).save(*args, **kwargs)

        # The md5 changes each time the object is saved
        #
        if not self.md5:
            key = '{}{}{}'.format(datetime.now(), self.pk, self.user)
            self.md5 = md5(key.encode('utf8')).hexdigest()

        super(TBUser, self).save(*args, **kwargs)


class TBAdminContact(TimeStampedModel):
    """
    People who receive system notification emails about file uploads, inquiries, etc
    """
    email = models.EmailField()
    name = models.CharField(max_length=255)

    is_active = models.BooleanField(default=True)

    def __str__(self):
        return self.email

    class Meta:
        verbose_name = 'TB Admin Contact'
        verbose_name_plural = 'TB Admin Contacts'
        ordering = ('email', )
