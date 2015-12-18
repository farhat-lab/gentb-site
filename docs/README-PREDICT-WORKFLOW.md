# Predict Workflow

There are multiple, somewhat disconnected steps in the Predict Workflow.

Overall, the Predict workflow steps are as follows, with notices regarding gaps:

  1. [Submission of Predict Form, including a Dropbox link](#1-submission-of-predict-form-including-a-dropbox-link)
  2. [Confirmation of Predict Information](#2-confirmation-of-predict-information)
  3. [Download Dropbox Files ("cron job")](#3-download-dropbox-files-cron-job)
  4. [Run Pipeline ("cron job")](#4-run-pipeline-cron-job)
  5. [Check for result files](#5-check-for-result-files)


## 1. Submission of Predict Form, including a Dropbox link

**Purpose:** Allow a registered user to submit a Dataset in the form of a Dropbox link as well as a description of the data.  Submission of a dropbox link containing VCF or FastQ files results in:
    - The creation of a PredictDataset object
    - The creation of a directory to hold files related to the PredictDataset

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

**Purpose:** Confirm that the PredictDataset information is accurate.  This page lists the files identified from the DropboxLink.

 - **Actions upon confirmation**
    - Change PredictDataset status to ```DATASET_STATUS_CONFIRMED```
 - **Actions upon canceling**
    - Delete PredictDataset and the associated file directory
    - File: [views_upload.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/predict/views_upload.py)

## 3. Download Dropbox Files ("cron job")

**Purpose:** Periodically (every 10 minutes), check if any new PredictDataset objects have been both added and confirmed.  If so, begin download the files in the appropriate directory for the PredictDataset.  (This directory name may be seen in the UI--the regular UI as well as Django admin)

  **To do**: If this process fails, the PredictDataset receives a status of [DATASET_STATUS_FILE_RETRIEVAL_ERROR](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/predict/models.py#L34).  However, an email should also be sent to the administrator.

  - This is activated re: a long running supervisord process -- e.g. it should be moved to a cron job once the cron server can use the needed libraries
    -  [run_dropbox_and_pipeline.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/cron_scripts/run_dropbox_and_pipeline.py)
      - For supervisord, search for ```[program:gentb_website]``` in the file [hms_supervisord.conf]( https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/tb_website/settings/hms_supervisord.conf)

  - **Additional dropbox info**
    - The download functionality requires the creation of a Dropbox app and associated access token: https://blogs.dropbox.com/developers/2015/08/new-api-endpoint-shared-link-metadata/
      - The token may be found in the [settings file](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/tb_website/settings/template_secret_settings.json) under ```DROPBOX_ACCESS_TOKEN```
        - On production, the settings file is named [```secret_settings_prod_hms.json```](https://github.com/IQSS/gentb-site/blob/master/docs/README-SETUP-ORCHESTRA.md#add-production-settings)
    - "Standalone" script that can be used to test the downloading of files from Dropbox links:
      - [dropbox_retriever.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/dropbox_helper/dropbox_retriever.py)
    - Wrapper object used to run the standalone script above:
      - [dropbox_retrieval_runner.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/dropbox_helper/dropbox_retrieval_runner.py)
    - This "wrapper":
      - Only tries downloads from [PredictDatast objects](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/predict/models.py)  with a status of [DATASET_STATUS_CONFIRMED](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/predict/models.py#L31)
      - Uses Dropbox links from [PredictDatast objects](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/predict/models.py) objects and updates the ```PredictDataset``` object statuses appropriately
      - Creates associated [DropboxRetrievalLog](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/dropbox_helper/models.py) objects--each one has an OneToOneField relationship with a  ```PredictDataset```


## 4. Run Pipeline ("cron job")

**Purpose:** Periodically (every 10 minutes), check if any new PredictDataset objects have finished downloading the related Dropbox files without error.  If so, run the pipeline script against the directory for this PredictDataset.

Preparation for this run includes the creation of a file named [gentb_status_feedback.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/templates/feedback/gentb_status_feedback.py) which is placed in the PredictDataset directory.  This python file contains callback information so that the pipeline script may send an end of job status message to the genTB server via API call.

  - Similar to the dropbox functionality, this is activated re: a long running supervisord process -- e.g. it should be moved to a cron job once the cron server can use the needed libraries
    -  [run_dropbox_and_pipeline.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/cron_scripts/run_dropbox_and_pipeline.py)
    - For supervisord, search for ```[program:gentb_website]``` in the file [hms_supervisord.conf]( https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/tb_website/settings/hms_supervisord.conf)

  - Searches for PredictDataset objects with a status of [DATASET_STATUS_FILE_RETRIEVAL_COMPLETE](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/predict/models.py#L35)
  - Kicks off the pipeline script using [pipeline_hardcoded_script_runner.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/predict/pipeline_hardcoded_script_runner.py)
    - Creates a [DatasetScriptRun object](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/predict/models.py) to record the results of the run--including the start time.
    - Constructs the appropriate pipeline command to run (VCF, FastQ single-end, FastQ pair-end, etc)
      - These are dynamically created depending on the data type and file extensions
  - This process is kicked off via a python Popen command--and then left to run without any further checking:
    - [script_runner_basic.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/script_helper/script_runner_basic.py)

## 5. Check for result files

**Purpose:** On a successful run, the pipeline scripts should create the output files described in [result_file_info.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/utils/result_file_info.py).  As of this writing (12/18/2015), these files include:
  - result.json   (predict script result)
  - matrix.csv    (pipeline analyze result)
  - heatmap.html  (predict script result)


The process is:
  1. As its last step the cluster script runs this python script [gentb_status_feedback.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/templates/feedback/gentb_status_feedback.py) which was created in the above step
  2. This script does the following:
    1. Makes a simple (somewhat naive) check to see if the results files (result.json, matrix.csv, etc):
       - (a) if these files exist and
       - (b) that they are more than 0 bytes in size.
    2. Make an API call to the genTB server to indicate success (true/false) and either an error message or the data contained in result.json
    3. These results are stored in [DatasetScriptRun objects](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/apps/predict/models.py)

Preparation for this run includes the creation of a file named [gentb_status_feedback.py](https://github.com/IQSS/gentb-site/blob/master/gentb_website/tb_website/templates/feedback/gentb_status_feedback.py) which is placed in the PredictDataset directory.  This python file contains callback information so that the pipeline script may send an end of job status message to the genTB server via API call.


fyi: This is a slightly dated diagram that captures the basic process--the numbering in the diagram __does not__ relate to the numbering in this document.
![predict workflow](images/predict-workflow.png?raw=true "Predict Workflow")
