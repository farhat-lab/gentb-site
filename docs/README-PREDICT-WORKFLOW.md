# Predict Workflow

There are multiple, somewhat disconnected steps in the Predict Workflow.

Overall, the Predict workflow steps are as follows, with notices regarding gaps:

  1. Submission of Predict Form, including a Dropbox link
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

## 3. Download Dropbox Files ("cron job")

## 4. Run Pipeline ("cron job")

## 5. Check for result files


![predict workflow](images/predict-workflow.png?raw=true "Predict Workflow")
