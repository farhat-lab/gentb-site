
from apps.predict.models import *

statuses = [
('DATASET_STATUS_NOT_READY_ID', 'Dataset NOT CONFIRMED for Dropbox file retrieval'),\
('DATASET_STATUS_CONFIRMED_ID', 'Dataset CONFIRMED for Dropbox file retrieval'),\
('DATASET_STATUS_FILE_RETRIEVAL_STARTED', 'Dropbox file retrieval STARTED'),\
('DATASET_STATUS_FILE_RETRIEVAL_ERROR', 'Dropbox file retrieval ERROR'),\
('DATASET_STATUS_FILE_RETRIEVAL_COMPLETE', 'Dropbox file retrieval COMPLETE'),\
('DATASET_STATUS_PROCESSING_STARTED_ID', 'Pipeline processing STARTED'),\
('DATASET_STATUS_PROCESSED_SUCCESS', 'Pipeline processing SUCCESS'),\
('DATASET_STATUS_PROCESSED_FAILED', 'Pipeline processing ERROR'),\
]

cnt = 0
for machine_name, human_name in statuses:
    cnt += 1
    (p, created) = PredictDatasetStatus.objects.get_or_create(id=cnt)#,\
    p.name = machine_name
    p.human_name = human_name
    p.sort_order = cnt
    #        name=name,\
    #        human_name=human_name,\
    #    sort_order=cnt)
    p.save()
