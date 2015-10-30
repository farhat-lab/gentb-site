

# Dropbox file retrieval.

## After confirm, new process kicked off to get file


1. command line script params
    - dropbox_url
    - callback_url
    - md5 of predict dataset, for callback
    - destination directory


(a) DropboxRetriever.py
(b) dropbox_retriever_runner.py -> command line args to wrap (a)
(c) dropbox_retrieve_util.py -> kick off (b)
