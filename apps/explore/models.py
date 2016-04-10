from django.db import models
from model_utils.models import TimeStampedModel
from django.utils.text import slugify


class ExploreDataFileInfo(TimeStampedModel):
    name = models.CharField(max_length=255, help_text='For internal use')
    active = models.BooleanField(default=True, help_text='The *most recently created* active entry will be used')

    # Link to download the Codebook
    codebook_file_url = models.URLField(help_text='Example: https://dataverse.harvard.edu/api/access/datafile/2694344')

    # Link to Two Ravens
    two_ravens_url = models.URLField(help_text='Example: https://rserve.dataverse.harvard.edu/dataexplore/gui.html?dfId=2693726&key=c54f07b7-5098-461c-adf3-a976c0d62f6e')

    def __str__(self):
        return '%s' % (self.name)

    class Meta:
        ordering = ('-created', )
        verbose_name = 'Explore Data File Information'
        verbose_name_plural = verbose_name

