"""
Submit the pipeline job as soon as the file download is complete.
"""
import sys
import time
from datetime import datetime

from django.core.management.base import BaseCommand
from django.contrib.sites.models import Site
from django.template.loader import get_template
from django.template import Context
from django.core.mail import send_mail
from django.conf import settings

from chore import JobSubmissionError

from apps.pipeline.models import PrepareError
from apps.predict.models import (
    PredictDataset, PredictStrain, get_timeout, STATUS_WAIT, STATUS_START, STATUS_ERROR
)
from apps.tb_users.models import TBAdminContact

# Which predict statuses should notify the user and admins?
STATUS_NOTIFY = (0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1)

def bitset(*args):
    """Turn a set of booleans into a bit array and then into an integer"""
    return int(''.join(reversed([str(int(bool(arg))) for arg in args])), 2)

def log(msg, *args, **kwargs):
    """Write consistantly to output"""
    kwargs['dt'] = datetime.now().isoformat()
    msg = "{dt}: " + msg
    sys.stderr.write(msg.format(*args, **kwargs) + "\n")

class Command(BaseCommand):
    """Schedule each of the pipeline tasks with the shell or LSF"""
    help = __doc__

    @staticmethod
    def submit_strain_pipeline(strain):
        """Submit a strain pipleine to the job queue"""
        status = strain.files_status
        if status in (STATUS_WAIT, STATUS_START):
            # Waiting or busy doing file upload
            return False
        if status is STATUS_ERROR:
            raise IOError("Download Error")

        try:
            return strain.run()
        except JobSubmissionError:
            raise ValueError("Can't run job")

    def notify_users(self):
        """Send emails to users if needed"""
        for predict in PredictDataset.objects.filter(has_notified=False):
            if STATUS_NOTIFY[predict.status]:
                self.send_email(predict.user, predict, 'predict/user_email.txt')
                for user in TBAdminContact.objects.filter(is_active=True):
                    if user != predict.user:
                        self.send_email(user, predict, 'predict/admin_email.txt')
                predict.has_notified = True
                predict.save()

    @staticmethod
    def send_email(user, predict, template_file):
        """Send a notification email to the given user about the predict"""
        if not user or not user.email:
            return

        subject = "GENTB PREDICT: {} [{}]".format(predict.title, predict.get_status())
        msg = get_template(template_file).render({
            'site': Site.objects.get_current(),
            'object': predict,
            'user': user,
        })
        send_mail(subject, msg, settings.DEFAULT_FROM_EMAIL, [user.email], fail_silently=False)

    def handle(self, **_):
        """Called from the command line"""
        # Limit all interactions to a timeout (usually a few weeks)
        qset = PredictStrain.objects.filter(dataset__created__gt=get_timeout())

        for strain in qset.filter(piperun__isnull=True):
            try:
                self.submit_strain_pipeline(strain)
                log("RUN: {} ({})", strain, strain.pipeline)
                time.sleep(0.25)
            except PrepareError as err:
                log("BAD: {}: {}".format(strain, err))
            except IOError as err:
                log("ERR: {} (Bad Download) {}", strain, err)
            except Exception as err: # pylint: disable=broad-except
                log("ERR: {} ({}): {}", strain, strain.pipeline, err)

        self.notify_users()
