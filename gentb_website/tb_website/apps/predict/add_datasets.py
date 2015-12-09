
from apps.predict.models import *

statuses = [
(False, 'DATASET_STATUS_NOT_READY', 'Dataset NOT CONFIRMED for Dropbox file retrieval'),\
(False, 'DATASET_STATUS_CONFIRMED', 'Dataset CONFIRMED for Dropbox file retrieval'),\
(False, 'DATASET_STATUS_FILE_RETRIEVAL_STARTED', 'Dropbox file retrieval STARTED'),\
(True, 'DATASET_STATUS_FILE_RETRIEVAL_ERROR', 'Dropbox file retrieval ERROR'),\
(False, 'DATASET_STATUS_FILE_RETRIEVAL_COMPLETE', 'Dropbox file retrieval COMPLETE'),\
(False, 'DATASET_STATUS_PROCESSING_STARTED', 'Pipeline processing STARTED'),\
(False, 'DATASET_STATUS_PROCESSED_SUCCESS', 'Pipeline processing SUCCESS'),\
(True, 'DATASET_STATUS_PROCESSED_FAILED', 'Pipeline processing ERROR'),\
]

cnt = 0
for is_error, machine_name, human_name in statuses:
    cnt += 1
    (p, created) = PredictDatasetStatus.objects.get_or_create(id=cnt)#,\
    p.is_error = is_error
    p.name = machine_name
    p.human_name = human_name
    p.sort_order = cnt
    #        name=name,\
    #        human_name=human_name,\
    #    sort_order=cnt)
    p.save()
