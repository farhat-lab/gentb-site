import os
from hashlib import md5
import json

from django.db import models
from django.core.files.storage import FileSystemStorage
from django.conf import settings
from datetime import datetime
from model_utils.models import TimeStampedModel

tb_file_system_storage = FileSystemStorage(location=settings.TB_SHARED_DATAFILE_DIRECTORY)

#def generate_new_filename(instance, filename):
#    #f, ext = os.path.splitext(filename)
#    instance.original_filename = basename(filename)
#    return join(instance.dataset.get_partial_path_for_datafile(), generate_storage_identifier())

class SharedFileInfo(TimeStampedModel):
    """
    Information from API call: https://api.github.com/repos/iqss/dataverse/milestones
    """
 
    first_name = models.CharField(max_length=75)
    middle_name =  models.CharField(max_length=20, blank=True)
    last_name = models.CharField(max_length=75)
    affiliation = models.CharField(max_length=255)
     
    contact_email = models.EmailField()
    
    title = models.CharField('Dataset title', max_length=255)

    description = models.TextField('Dataset description')

    file_obj = models.FileField('Data upload', upload_to='shared-files/%Y/%m'\
                            , storage=tb_file_system_storage)

    has_prediction = models.BooleanField('Has prediction results?',default=False, help_text='auto-filled on save')
    prediction_results = models.TextField( blank=True)

    md5 = models.CharField(max_length=40, blank=True, db_index=True, help_text='auto-filled on save')

    def __unicode__(self):
        return '%s' % (self.title)
    
    def __str__(self):
        return self.__unicode__()


    def set_prediction_results(self, results_dict):
        """
        Converts a python dict to JSON and sets the "prediction_results" attribute
        Note, this does not save the object
        """
        if not (type(results_dict) in [tuple, dict, list]):
            #return False
            raise TypeError('Expected a python dict, list, or tuple')
    
        try:
            json_str = json.dumps(results_dict)
        except:
            #return False
            raise Exception('Failed to convert python object to json. (type: "%s")' % type(results_dict))

        self.prediction_results  = json_str
        return True

    def get_prediction_results(self):
        if not self.prediction_results:
            return ''

        try:
            return json.loads(self.prediction_results)   # JSON string -> python; e.g. '[1, 2, 3]' -> [1, 2, 3]
        except:
            return ''
            #raise Exception('Failed to convert prediction result json to dict')



    def save(self, *args, **kwargs):
        if not self.id:
            super(SharedFileInfo, self).save(*args, **kwargs)

        self.md5 = md5('%s%s%s' % (str(datetime.now()), self.id, self.title)).hexdigest()

        if self.prediction_results:
            self.has_prediction = True
        else:
            self.has_prediction = False

        super(SharedFileInfo, self).save(*args, **kwargs)
        
    def filename(self):
        return os.path.basename(self.file_obj.name)
    
    class Meta:
        ordering = ('title', '-created')
        verbose_name_plural = 'Shared File Information'