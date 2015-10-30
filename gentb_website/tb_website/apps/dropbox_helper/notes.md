

# Dropbox file retrieval.

## After confirm, new process kicked off to get file


1. command line script params
    - dropbox_url
    - callback_url
    - md5 of predict dataset, for callback
    - destination directory


(A) dropbox_retriever - DropboxRetriever
    - Main class to retrieve files from a dropbox link

(B) dropbox_retriever_runner.py
    - Accepts command line arguments, including callback_url to run the DropboxRetriever
    - Uses callback_url to send results back to DropboxRetrievalLog

(C) dropbox_retrieve_util.py -> kick off (b)
    - Business logic to:
        1. Update status of the predict_dataset
        2. Kick of (b) dropbox_retriever_runner with command line args
