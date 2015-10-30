from datetime import datetime
import json
import os, sys, shutil
from os.path import basename, dirname, join, isdir, isfile, realpath
import re
import requests
import urllib2
import zipfile

if __name__=='__main__':
    # For local testing....
    django_dir = dirname(dirname(dirname(realpath(__file__))))
    sys.path.append(django_dir)
    os.environ['DJANGO_SETTINGS_MODULE'] = 'tb_website.settings.local'

from django.conf import settings

GENTB_FILE_PATTERNS = ['\.fastq$', '\.fastq\.', '\.vcf$', '\.vcf\.', ]
#GENTB_FILE_PATTERNS = ['\.fastq$', '\.fastq\.', '\.vcf$', '\.vcf\.', '\.txt$']

class DropboxRetriever:
    """
    Given a Dropbox shared link, scan it for matching files and retrieve them.

    This functionality uses the Dropbox Core API to retrieve metadata from a shared link.
        https://blogs.dropbox.com/developers/2015/08/new-api-endpoint-shared-link-metadata/
    This requires a Dataverse app (either app_key + secret OR or access_token)
    """
    def __init__(self, dbox_link, destination_dir, file_patterns=GENTB_FILE_PATTERNS):

        self.dbox_link = dbox_link
        self.is_directory_link = False

        self.destination_dir = destination_dir
        self.file_patterns = file_patterns

        # For Step 1 - retrieve metadata
        self.dropbox_link_metadata = None

        # For Step 2 - check file matches
        self.has_file_matches = False
        self.matching_files_metadata = []

        # For Step 3 - retrieve files
        self.target_zip_fullname = None
        self.final_file_paths = []

        # Error holders
        self.err_found = False
        self.err_msg = None

        self.initial_check()

    def add_err_msg(self, m):
        self.err_found = True
        self.err_msg = m

    def initial_check(self):

        if not self.is_valid_dropbox_link():
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

    def is_valid_dropbox_link(self):
        """
        Starts with "https://www.dropbox.com"
        Ends with: "?dl=0" or "?dl=1"

        Always changes ?dl=0" -> "?dl=1"
        """
        if self.dbox_link is None:
            self.add_err_msg("dbox_link cannot be None")
            return False

        d = self.dbox_link.lower()

        dbox_start = 'https://www.dropbox.com'
        if not d.startswith(dbox_start):
            self.add_err("Not a dropbox url.  Doesn't start with '%s'" % dbox_start)
            return False

        if d.endswith('?dl=0'):
            self.dbox_link = self.dbox_link.replace('?dl=0', '?dl=1')
            return True

        if d.endswith('?dl=1'):
            return True

        self.add_err('Does not appear to be a valid download link.  Should end with "?dl=0" or "?dl=1"')
        return False


    def step1_retrieve_metadata(self):
        """
        Retrieve metadata about this dropbox link
        """
        # Prepare request params (the dropbox link)
        #
        params = dict(link=self.dbox_link)

        # Add token to the header
        #
        headers = {'Authorization': 'Bearer {0}'.format(settings.DROPBOX_ACCESS_TOKEN)}

        # Make the request
        #
        r = requests.post('https://api.dropbox.com/1/metadata/link',
                        data=params,
                        headers=headers
                        )
        if r.status_code != 200:
            dbox_err_msg = None
            try:
                rjson = r.json()
                if 'error' in rjson:
                    emsg = 'Error: {0}'.format(rjson.get('error'))
                    if r.status_code == 403:
                        dbox_err_msg = 'The dropbox link did not work.  Please check the url. ({0})'.format(emsg)
                    else:
                        dbox_err_msg = 'The dropbox link returned an error. Status code: {0}\n{1}'.format(r.status_code, emsg)
            except:
                dbox_err_msg = 'The dropbox link returned an error. Status code: {0}\n{1}'.format(r.status_code, emsg)

            self.add_err_msg(dbox_err_msg)
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

        # -------------------------------------
        # Handle a Dropbox link to a file?
        # -------------------------------------
        if self.dropbox_link_metadata['is_dir'] is False:   # Yes, this is a file
            self.is_directory_link = False

            # Does it match?
            fpath = self.dropbox_link_metadata['path']
            if self.does_file_match_criteria(fpath):
                self.matching_files_metadata.append(fpath)
                return True
            else:
                self.add_err_msg('No files match what we are looking for.\
                        Please make sure you have at least one ".fastq" or ".vcf" file.')
                return False

        # -------------------------------------
        # Handle a Dropbox link to a directory?
        # -------------------------------------
        self.is_directory_link = True

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
            self.add_err_msg('No files match what we are looking for. Please make sure you have at least one ".fastq" or ".vcf" file.')
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
        Download the dropbox link as a .zip.
        There is currently no method download selected files from a shared dierctory link.
        Desired paths stored in "self.matching_files_metadata"
        """
        if self.err_found:
            return False

        assert isinstance(self.matching_files_metadata, list), \
            "Do not call this unless dropbox metadata (step 2) has files that we want"
        assert len(self.matching_files_metadata) > 0, \
            "No matching files. Do not call this unless dropbox metadata (step 2) has files that we want"

        if self.is_directory_link:
            if self.step4a_download_directory_as_zip():
                if self.step4b_organize_files():
                    return True
        else:
            self.step3a_download_single_file()
            return True
        return False

    def step3a_download_single_file(self):

        if self.err_found:
            return False

        # -------------------------------------
        # Set the .zip name and destination
        # -------------------------------------
        target_fname = 'dropbox_download_{0}.zip'.format(datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))

        self.target_fullname = join(self.destination_dir, basename(self.dropbox_link_metadata['path'][1:]))

        chunk_size = 16 * 1024
        req = urllib2.urlopen(self.dbox_link)

        # -------------------------------------
        # Download it!
        # -------------------------------------
        print 'pre download: {0} @ {1}'.format(self.dbox_link, datetime.now())
        with open(self.target_fullname, 'wb') as fp:
            print 'downloading...'
            shutil.copyfileobj(req, fp, chunk_size)

        print 'file written: {0} @ {1}'.format(self.target_fullname, datetime.now())
        return True


    def step4a_download_directory_as_zip(self):
        if self.err_found:
            return False

        # -------------------------------------
        # Set the .zip name and destination
        # -------------------------------------
        target_fname = 'dropbox_download_{0}.zip'.format(datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
        self.target_zip_fullname = join(self.destination_dir, target_fname)

        chunk_size = 16 * 1024
        req = urllib2.urlopen(self.dbox_link)

        # -------------------------------------
        # Download it!
        # -------------------------------------
        print 'pre download: {0} @ {1}'.format(self.dbox_link, datetime.now())
        with open(self.target_zip_fullname, 'wb') as fp:
            print 'downloading...'
            shutil.copyfileobj(req, fp, chunk_size)

        print 'file written: {0} @ {1}'.format(self.target_zip_fullname, datetime.now())

        print '-' * 50
        for fp in self.final_file_paths:
            print fp

        return True


    def step4b_organize_files(self):

        if not isfile(self.target_zip_fullname):
            self.add_err_msg("The dropbox .zip file was not found: {0}".format(self.target_zip_fullname))
            return False

        # -------------------------------------
        # Unzip the download
        # -------------------------------------
        unzip_dir = join(self.destination_dir, self.dropbox_link_metadata['path'][1:])
        try:
            dbox_zip = zipfile.ZipFile(self.target_zip_fullname, 'r')
            dbox_zip.extractall(unzip_dir)
            dbox_zip.close()
        except:
            self.add_err_msg("Failed to unzip file '{0}' to directory '{1}'".format(\
                                self.target_zip_fullname, self.destination_dir))
            return False

        # -------------------------------------
        # Remove directories and unwanted files
        # -------------------------------------
        for item in os.listdir(unzip_dir):
            full_item = join(unzip_dir, item)

            # -------------------------------------
            # The files must be at the first level
            # delete any embedded directories
            # -------------------------------------
            if isdir(full_item):
                shutil.rmtree(full_item)
                print 'directory removed: {0}'.format(full_item)
                continue    # go to the next item

            # -------------------------------------
            # It's a file, one of ours?
            # -------------------------------------
            needed_file = False
            for fname in dr.matching_files_metadata:
                if full_item.endswith(fname):
                    needed_file = True
            if not needed_file:
                os.remove(full_item)
                print 'file removed: {0}'.format(full_item)
            else:
                self.final_file_paths.append(full_item)
                print 'file kept: {0}'.format(full_item)

        if len(self.final_file_paths) == 0:
            self.add_err_msg("The desired files were not downloaded")
            return False

        print '-' * 50
        for fp in self.final_file_paths:
            print fp

        return True

if __name__=='__main__':
    # For local testing....see sys.path at top if using this non-locally

    #example_dlink = 'https://www.dropbox.com/s/4tqczonkaeakvua/001.txt?dl=0'
    #example_dlink = 'https://www.dropbox.com/sh/vbicdol2e8mn57r/AACeJzBhUpgxTNjj6jHL2UJoa?dl=0'
    # dir of files
    example_dlink = 'https://www.dropbox.com/sh/19krhpbo4ph93rp/AAB6z3SpKs3w7jHy0bVi4JtPa?dl=0'
    # single file
    example_dlink = 'https://www.dropbox.com/s/tipb48hahtmbk05/008.1.fastq.txt?dl=0'
    dest_dir = '/Users/rmp553/Documents/iqss-git/gentb-site/scratch-work/test-files'

    # Initialize
    #
    dr = DropboxRetriever(example_dlink, dest_dir, GENTB_FILE_PATTERNS)
    if dr.err_found:
        print dr.err_msg
        sys.exit(1)

    # Get the metadata
    #
    if not dr.step1_retrieve_metadata():
        print dr.err_msg
        sys.exit(1)

    # Does it have what we want?
    #
    if not dr.step2_check_file_matches():
        print dr.err_msg
        sys.exit(1)

    print dr.matching_files_metadata
    #sys.exit(1)

    # Download the files
    #
    if not dr.step3_retrieve_files():
        print dr.err_msg
        sys.exit(1)
