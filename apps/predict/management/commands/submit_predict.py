"""
Submit the pipeline job as soon as the file download is complete.
"""
import os
import sys
import time
import shutil
from datetime import datetime

from django.core.management.base import BaseCommand
from django.contrib.sites.models import Site
from django.template.loader import get_template
from django.template import Context
from django.core.mail import send_mail
from django.conf import settings

from chore import JobSubmissionError

from apps.pipeline.models import PrepareError
from apps.predict.models import PredictDataset, PredictStrain, get_timeout
from apps.tb_users.models import TBAdminContact

# Which predict statuses should notify the user and admins?
STATUS_NOTIFY = (
    'FILE_ERROR', 'RUN_ERROR', 'RUN_DONE', 'READY', 'INVALID', 'TIMEOUT',
)

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
        if status in ('FILE_START', 'FILE_WAIT', 'FILE_NONE'):
            # Waiting or busy doing file upload
            return False
        if status is 'FILE_ERROR':
            raise IOError("Download Error")

        for st_file in (strain.file_one, strain.file_two):
            if not st_file:
                continue
            if st_file.filename.endswith('vcf.gz'):
                st_file.decompress()
            elif st_file.filename.endswith('.fq'):
                st_file.rename(st_file.filename[:-3] + ".fastq")
            elif st_file.filename.endswith('fq.gz'):
                st_file.rename(st_file.filename[:-3] + ".fastq.gz")

        try:
            if not strain.run():
                raise ValueError("Run didn't work")
        except JobSubmissionError:
            raise ValueError("Can't run job")

    def notify_users(self):
        """Send emails to users if needed"""
        for predict in PredictDataset.objects.filter(has_notified=False):
            if predict.status in STATUS_NOTIFY:
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

        for strain in qset.filter(piperun__isnull=True, pipeline__isnull=False, pipeline__disabled=False):
            try:
                if not self.submit_strain_pipeline(strain):
                    continue
                log("RUN: {} ({})", strain, strain.pipeline)
                time.sleep(0.25)
            except PrepareError as err:
                log("BAD: {}: {}".format(strain, err))
                pipeline = strain.pipeline
                pipeline.disabled = True
                pipeline.save()
            except IOError as err:
                log("ERR: {} (Bad Download) {}", strain, err)
            except Exception as err: # pylint: disable=broad-except
                log("ERR: {} ({}): {}", strain, strain.pipeline, err)

        clean_predict_dir()
        self.notify_users()


def clean_predict_dir():
    """
    Attempt to remove any of the directories no longer in use to save space.
    """
    parent = '/n/groups/gentb_www/predictData'
    if not os.path.isdir(parent):
        return

    log("Starting cleaning process: {}".format(parent))
    expected = list(PredictDataset.objects.values_list('file_directory', flat=True))
    found = os.listdir(parent)

    for fname in found:
        path = str(os.path.join(parent, fname))
        if path not in expected:
            try:
                shutil.rmtree(path)
            except OSError:
                sys.stderr.write(" [!] {}\n".format(path))
            else:
                sys.stderr.write(" [X] {}\n".format(path))

    for path in expected:
        if not os.path.isdir(path):
            sys.stderr.write(" [R] {}\n".format(path))
            PredictDataset.objects.filter(file_directory=path).delete()
