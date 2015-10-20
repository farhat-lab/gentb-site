import json
from os.path import isdir
from pprint import pprint
import re
import requests
import sys
from dropbox_info import DROPBOX_ACCESS_TOKEN

#GENTB_FILE_PATTERNS = ['\.fastq$', '\.fastq\.', '\.vcf$', '\.vcf\.', ]
GENTB_FILE_PATTERNS = ['\.fastq$', '\.fastq\.', '\.vcf$', '\.vcf\.', '\.txt$']

class DropboxRetriever:
    """
    Given a Dropbox shared link, scan it for matching files and retrieve them.

    This functionality uses the Dropbox Core API to retrieve metadata from a shared link.
        https://blogs.dropbox.com/developers/2015/08/new-api-endpoint-shared-link-metadata/
    This requires a Dataverse app (either app_key + secret OR or access_token)
    """
    def __init__(self, dbox_link, destination_dir, file_patterns=[]):

        self.dbox_link = dbox_link
        self.destination_dir = destination_dir
        self.file_patterns = file_patterns

        # For Step 1 - retrieve metadata
        self.dropbox_link_metadata = None

        # For Step 2 - check file matches
        self.has_file_matches = False
        self.matching_files_metadata = []

        # For Step 3 - retrieve files
        self.file_paths = []

        # Error holders
        self.err_found = False
        self.err_msg = None

        self.initial_check()

    def add_err_msg(self, m):
        self.err_found = True
        self.err_msg = m

    def initial_check(self):

        if self.dbox_link is None:
            self.add_err_msg("dbox_link cannot be None")
            return False

        if self.destination_dir is None:
            self.add_err_msg("destination_dir cannot be None")
            return False

        if not isdir(self.destination_dir):
            self.add_err_msg("The destination directory does not exist: {0}".format(self.destination_dir))
            return False

        if self.file_patterns is None:
            self.file_patterns = []

        return True

    def step1_retrieve_metadata(self):
        """
        Retrieve metadata about this dropbox link
        """
        global DROPBOX_ACCESS_TOKEN

        # Prepare request params (the dropbox link)
        #
        params = dict(link=self.dbox_link)

        # Add token to the header
        #
        headers = {'Authorization': 'Bearer {0}'.format(DROPBOX_ACCESS_TOKEN)}

        # Make the request
        #
        r = requests.post('https://api.dropbox.com/1/metadata/link',
                        data=params,
                        headers=headers
                        )
        if r.status_code != 200:
            self.add_err_msg("The dropbox request returned an error. Status code: {0}\nError: {1}".format(r.status_code, r.text))
            return False

        try:
            rjson = r.json()
        except:
            self.add_err_msg("Failed to turn link metadata from dropbox to JSON: {0}".format(r.text))
            return False

        # Set object variable
        #
        self.dropbox_link_metadata = rjson

        # Does this look correct?
        if not "is_dir" in self.dropbox_link_metadata:
            self.add_err_msg("Metadata doesn't look good. Missing key 'is_dir': {0}".format(r.text))
            return False

        #pprint(self.dropbox_link_metadata, depth=4)
        print json.dumps(self.dropbox_link_metadata, indent=4, sort_keys=True)
        return True


    def step2_check_file_matches(self):
        """
        Does the dropbox metadata contain files that we're looking for?
        """
        if self.err_found:
            return False
        assert isinstance(self.dropbox_link_metadata, dict),\
            "Do not call this unless step 1, retrieval of Dropbox metadata is successful"
        assert "is_dir" in self.dropbox_link_metadata,\
            "Do not call this unless step 1, retrieval of Dropbox metadata is successful. (Metadata doesn't look good. Missing key 'is_dir')"

        # Is this link to a file?
        if self.dropbox_link_metadata['is_dir'] is False:   # Yes, this is a file
            # Does it match?
            fpath = self.dropbox_link_metadata['path']
            if self.does_file_match_criteria(fpath):
                self.matching_files_metadata.append(fpath)
                return True
            else:
                self.add_err_msg('No files match what we are looking for.')
                return False

        # This is a directory, check each file
        #
        if not 'contents' in self.dropbox_link_metadata:
            self.add_err_msg("The directory is empty.  No 'contents' key")
            return False

        # Iterate through files
        #
        for file_info in self.dropbox_link_metadata['contents']:
            if file_info['is_dir'] is True: # This is a subdirectory, keep going
                continue
            fpath = file_info['path']   # Assume the 'path' is always here
            if self.does_file_match_criteria(fpath):
                self.matching_files_metadata.append(fpath)

        if len(self.matching_files_metadata) == 0:
            self.add_err_msg('No files match what we are looking for.')
            return False

        # We've got something
        return True


    def does_file_match_criteria(self, fpath):
        """
        For genTB, check if this has a .fastq or .vcf extension
        """
        if fpath is None:
            return False

        if len(self.file_patterns) == 0:    # nothing to check
            return True

        # Iterate through the file patterns and see if one matches
        #
        for pat in self.file_patterns:
            if re.search(pat, fpath, re.IGNORECASE):
                return True
        return False


    def step3_retrieve_files(self):
        """
        Download the dropbox file paths stored in "self.matching_files_metadata"
        """
        if self.err_found:
            return False
        assert isinstance(self.matching_files_metadata, []), \
            "Do not call this unless dropbox metadata (step 2) has files that we want"
        assert len(self.matching_files_metadata) > 0, \
            "Do not call this unless dropbox metadata (step 2) has files that we want"

        

print type(self.dropbox_link_metadata)
assert isinstance(self.dropbox_link_metadata, dict),\
    "Do not call this unless step 1, retrieval of Dropbox metadata is successful"

        if self.err_found:
            return False

        pass

if __name__=='__main__':
    #example_dlink = 'https://www.dropbox.com/s/4tqczonkaeakvua/001.txt?dl=0'
    example_dlink = 'https://www.dropbox.com/sh/vbicdol2e8mn57r/AACeJzBhUpgxTNjj6jHL2UJoa?dl=0'
    dest_dir = '/Users/rmp553/Documents/iqss-git/gentb-site/scratch-work/test-files'

    # initiate
    #
    dr = DropboxRetriever(example_dlink, dest_dir, GENTB_FILE_PATTERNS)
    if dr.err_found:
        print dr.err_msg
        sys.exit(1)

    if not dr.step1_retrieve_metadata():
        print dr.err_msg
        sys.exit(1)

    if not dr.step2_check_file_matches():
        print dr.err_msg
        sys.exit(1)

    print dr.matching_files_metadata
