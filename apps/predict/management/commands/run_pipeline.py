
# pylint: disable=no-init,no-self-use
"""
Run the pipeline.
"""

from django.core.management.base import BaseCommand, CommandError

from apps.predict.models import LOGGER, PredictDataset

class Command(BaseCommand):
    """Fires off the perl script to kick off the pipeline"""
    help = """
This fires off the perl script to kick off the pipeline for either
FastQ or VCF file analysis.
"""

    def add_arguments(self, parser):
        """No arguments for this command yet"""
        pass

    """
    Given a PredictDataset, run the appropriate
    cluster pipeline script to process it.
    """
    err_found = False
    err_message_title = None
    err_message = None

    def handle(self, **options):
        """
        (1) Retrieve the first PredictDataset with a status of:
            DATASET_STATUS_FILE_RETRIEVAL_COMPLETE
        (2) If such a dataset exists, run it through the pipeline
        """
        LOGGER.info("Run pipeline check: next dataset")

        # get some Dataset
        dataset = PredictDataset.objects.filter(\
                status=PredictDataset.STATUS_FILE_RETRIEVAL_COMPLETE).first()

        if dataset is None:
            return LOGGER.info("Nothing to check")

        # Run script
        LOGGER.info("Run pipeline for dataset: %s (%s)", dataset, dataset.pk)

        (ret, msg) = dataset.run_command()
        if not ret:
            LOGGER.error(str(msg))
        LOGGER.info("OK")


