

# Dropbox file retrieval.

## After confirm, cron job will catch new links and retrieve files


(A) dropbox_retriever - DropboxRetriever
    - Main class to retrieve files from a dropbox link

(B) dropbox_retriever_runner.py
    - Kicked off with a cron job
    - Checks recently upload dropbox links and tries to retrive their files
