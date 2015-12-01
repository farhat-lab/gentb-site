"""
This script checks the local "output" directory and sends feedback to the server
about the job.

/output/
    - result.json
    - matrix.csv
"""
from os.path import dirname, isdir, isfile, getsize, realpath
import os, sys

#class GenTBStatusFeedback:



#os.path.dirname(os.path.realpath(__file__))
