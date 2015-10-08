from django.conf import settings
from django.core.mail import send_mail

from apps.tb_users.models import TBAdminContact

def get_tb_admin_emails():
    """
    Pull email values from TB_ADMINS

    settings.TB_ADMINS = (
        ('Your Name', 'your_email@example.com'),
    )
    :return:
    """
    return [ info[1] for info in settings.TB_ADMINS ]

def send_mail_to_user_and_admins(subject, user_email, text_msg, html_message=None):
    assert user_email is not None, 'user_email cannot be None'
    assert subject is not None, 'subject cannot be None'
    assert text_msg is not None, 'text_msg cannot be None'

    # (1) check the database
    to_emails = TBAdminContact.objects.filter(is_active=True).values_list('email', flat=True)

    # (2) default to the TB_ADMINS in the settings file
    if to_emails.count() == 0:
        to_emails = get_tb_admin_emails()
    else:
        to_emails = list(to_emails)

    if to_emails is None or len(to_emails) == 0:
        raise ValueError('settings.TB_ADMINS must have at least one contact email address.  Check settings/local.py or settings/production.py')

    to_emails.append(user_email)

    send_mail(subject,
              text_msg,
              settings.DEFAULT_FROM_EMAIL,
              to_emails,
              fail_silently=False,
              html_message=html_message)

def send_email_to_admins(subject, text_msg, html_message=None):
    assert subject is not None, 'subject cannot be None'
    assert text_msg is not None, 'text_msg cannot be None'

    # (1) check the database
    to_emails = TBAdminContact.objects.filter(is_active=True).values_list('email', flat=True)

    # (2) default to the TB_ADMINS in the settings file
    if to_emails.count() == 0:
        to_emails = get_tb_admin_emails()

    if to_emails is None or len(to_emails) == 0:
        raise ValueError('settings.TB_ADMINS must have at least one contact email address.  Check settings/local.py or settings/production.py')

    send_mail(subject,
              text_msg,
              settings.DEFAULT_FROM_EMAIL,
              to_emails,
              fail_silently=False,
              html_message=html_message)

