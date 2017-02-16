from django.template.loader import render_to_string
from apps.utils.email_util import send_email_to_admins, send_mail_to_user_and_admins

from .utils import get_site_url

def send_new_dataset_message_to_tb_admins(dataset):
    """
    After a prediction file is uploaded, send the TB ADMINS a message

    :param dataset:
    :return:
    """
    subject = "genTB: New file uploaded"

    d = dict(dataset=dataset,
             tb_user=dataset.user,
             subject=subject,
             SITE_URL=get_site_url())

    text_message = render_to_string('predict/email/notify_file_upload.txt', d)
    html_message = render_to_string('predict/email/notify_file_upload.html', d)

    send_email_to_admins(subject, text_message, html_message)



def send_dataset_run_message_to_tb_admins_and_user(dataset_script_run):

    if dataset_script_run.result_success:
        subject = "genTB: New file processed (Success)"
    else:
        subject = "genTB: New file processed (Fail)"

    user_email = dataset_script_run.dataset.user.email

    d = dict(dataset=dataset_script_run.dataset,
             tb_user=dataset_script_run.dataset.user,
             script_run=dataset_script_run,
             subject=subject,
             SITE_URL=get_site_url())

    if dataset_script_run.result_success:
        text_message = render_to_string('predict/email/pipeline_success_run.txt', d)
        html_message = render_to_string('predict/email/pipeline_success_run.html', d)
    else:
        text_message = render_to_string('predict/email/pipeline_fail_run.txt', d)
        html_message = render_to_string('predict/email/pipeline_fail_run.html', d)

    send_mail_to_user_and_admins(subject, user_email, text_message, html_message)

