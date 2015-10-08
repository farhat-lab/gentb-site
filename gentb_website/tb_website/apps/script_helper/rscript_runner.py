from __future__ import print_function
import os
import sys, json
import subprocess # for running rscript
from apps.utils.msg_util import *

from django.conf import settings  

class RScriptRunner:
    
    def __init__(self, file_fullpath):
        self.CURRENT_DIR = os.getcwd()
        
        assert file_fullpath, not None
        
        #msgt('file to check: %s' % file_fullpath)
        
        self.rcommand_base = 'Rscript TBPredict.R'
        self.file_fullpath = file_fullpath
        self.err_found = False
        self.err_msg_list = []
        self.formatted_response = None      # python dict if successful
        
    def add_err_msg(self, err_msg):
        self.err_found = True
        if err_msg:
            self.err_msg_list.append(err_msg)
    
    def get_err_msgs(self, as_html=True):
        assert type(self.err_msg_list), list
        
        if len(self.err_msg_list) == 0:
            return None
            
        if as_html:
            delim = '<br />'
        else:
            delim = '\n'

        return delim.join(self.err_msg_list)
    
    def return_to_previous_working_directory(self):
        os.chdir(self.CURRENT_DIR)      # return to previous working directory
        
        
    def run_command_on_file(self):
        """
        Run the R script and do lots of error checking.

        If the script runs successfully, the result is stored as a python dict in self.formatted_response
        """
        # Does the file exist
        if not os.path.isfile(self.file_fullpath):
            self.add_err_msg('The file was not found: "%s"' % os.path.basename(self.file_fullpath))
            return False
            
        # Format R command
        rscript_full_command = '%s %s' % (self.rcommand_base, self.file_fullpath)
        #rscript_results = subprocess.check_output(rscript_full_command.split())
        
        # Go to the R scripts directory
        #
        os.chdir(settings.R_SCRIPTS_PATH)
        msg('change os.chdir: %s' % settings.R_SCRIPTS_PATH)
        msgt('rscript_full_command: %s' % rscript_full_command)
        try:
            rscript_process = subprocess.Popen(rscript_full_command.split(),
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE
                                    )
            stdout_value, stderr_value = rscript_process.communicate() #proc.communicate
            rscript_results = stdout_value
            rscript_err_results = stderr_value

            #rscript_results = rscript_process.communicate()[0]
            #rscript_results = subprocess.check_output(rscript_full_command.split())
            msg ('rscript_results: %s' % rscript_results)
            msg ('rscript_err_results: %s' % rscript_err_results)

        except OSError as e:
            self.add_err_msg('Sorry!  There was was an error running the script. (1)')
            self.add_err_msg('OSError: %s' % e)
            self.return_to_previous_working_directory()

            return False
        except subprocess.CalledProcessError as e:

            self.add_err_msg('Sorry!  There was was an error running the script. (2)')
            self.add_err_msg("CalledProcess error: %s" % e.message)

            self.return_to_previous_working_directory()
            return False
        
        except:
            e = sys.exc_info()[0]
            self.add_err_msg('Sorry!  There was an error running the script. (3)')
            if len(str(e)) > 0:
                self.add_err_msg('Error: %s' % e)
            self.add_err_msg("(file name: %s)" % os.path.basename(self.file_fullpath))
            self.return_to_previous_working_directory()
            return False
        
        self.return_to_previous_working_directory()

        if rscript_err_results and rscript_err_results.find('Fatal error:') > -1:
            self.add_err_msg('Sorry!  There was an error running the script. (4)')
            self.add_err_msg("Error: %s" % rscript_err_results)
            return False

        if rscript_results and rscript_results.find('Fatal error:') > -1:
            self.add_err_msg('Sorry!  There was an error running the script. (5)')
            self.add_err_msg("Error: %s" % rscript_results)
            self.add_err_msg("R script path: %s" % settings.R_SCRIPTS_PATH)

            return False

        json_resp = None
        try:
            rscript_results = rscript_results.strip()
            json_resp = json.loads(rscript_results)
        except ValueError as e:
            self.add_err_msg('Sorry!  There was an error processing the script results. (6)')
            self.add_err_msg("\n(ValueError: %s)" % e)
            self.add_err_msg("\n(script results: %s)" % rscript_results)
            #self.add_err_msg("Error({0}): {1}".format(e.errno, e.strerror))
            return False        
        except:
            e = sys.exc_info()[0]
            self.add_err_msg('Sorry!  There was an error running the script. (7)')
            self.add_err_msg('Error: %s' % e)
            self.add_err_msg("(file name: %s)" % os.path.basename(self.file_fullpath))
            return False
        
        self.formatted_response = json_resp
        return True


    def get_formatted_response(self):
        if self.err_found:
            return None
        return self.formatted_response
        