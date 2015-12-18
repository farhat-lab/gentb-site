# Predict Workflow

There are multiple, somewhat disconnected steps in the Predict Workflow.

Overall, the Predict workflow steps are as follows, with notices regarding gaps:

  1. [Submission of Predict Form, including a Dropbox link](#1-submission-of-predict-form-including-a-dropbox-link)
  2. Confirmation of Predict Information
  3. Download Dropbox Files ("cron job")
  4. Run Pipeline ("cron job")
  5. Check for result files


## 1. Submission of Predict Form, including a Dropbox link

- **Error checking notes**
  - The Dropbox link Metadata is checked to see if it contains
    the chosen file types, VCF or FastQ. This is part of the Django forms's
    error check.
    - Form name: [UploadPredictionDataForm](../gentb_website/tb_website/apps/predict/forms.py#L47)
       - Method: clean() method, see ```get_dropbox_metadata_from_link```
         - from: ```from apps.dropbox_helper.dropbox_util import get_dropbox_metadata_from_link```
- **Actions upon success**
  - Creation of a new ```PredictDataset object``` and related ```DropboxRetrievalLog``` object
    - Includes the creation of a file directory for the ```PredictDataset object```
        - See [PredictDataset](../gentb_website/tb_website/apps/predict/models.py)
        - Method: ```create_dataset_directory_name```
  - See [UploadPredictionDataForm](../gentb_website/tb_website/apps/predict/forms.py)
    - Method: get_dataset(tb_user)

## 2. Confirmation of Predict Information
 - **Actions upon confirmation**
    - Change PredictDataset status to ```DATASET_STATUS_CONFIRMED```
 - **Actions upon canceling**
    - Delete PredictDataset and the associated file directory
    - File: [views_upload.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/predict/views_upload.py)

## 3. Download Dropbox Files ("cron job")
  - This is activated re: a long running supervisord process -- e.g. it should be moved to a cron job once the cron server can use the needed libraries
    -  [run_dropbox_and_pipeline.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/cron_scripts/run_dropbox_and_pipeline.py)
      - For supervisord, search for ```[program:gentb_website]``` [hms_supervisord.conf]( https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/tb_website/settings/hms_supervisord.conf)

  - **Additional dropbox info***    
    - The download functionality requires the creation of a Dropbox app and associated access token: https://blogs.dropbox.com/developers/2015/08/new-api-endpoint-shared-link-metadata/
      - The token may be found in the [settings file](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/tb_website/settings/template_secret_settings.json) under ```DROPBOX_ACCESS_TOKEN```
        - On production, the settings file is named [```secret_settings_prod_hms.json```](https://github.com/IQSS/gentb-site/blob/master/docs/README-SETUP-ORCHESTRA.md#add-production-settings)
    - "Standalone" script that can be used to test the downloading of files from Dropbox links:
      - [dropbox_retriever.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/dropbox_helper/dropbox_retriever.py)
    - Wrapper for the standalone script above.  This "wrapper":
      - Only tries downloads from PredictDatast objects with a status of [DATASET_STATUS_CONFIRMED](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/predict/models.py#L31)
      - Uses Dropbox links from ```PredictDataset``` objects and updates the ```PredictDataset``` object statuses appropriately
      - Creates associated ```DropboxRetrievalLog``` objects--each one has an OneToOneField relationship with a  ```PredictDataset```


## 4. Run Pipeline ("cron job")

## 5. Check for result files


![predict workflow](images/predict-workflow.png?raw=true "Predict Workflow")
