
from os.path import isfile, join, dirname, realpath
import os, sys

CURRENT_DIR = dirname(realpath(__file__))
RDIR = realpath(join(CURRENT_DIR, '../../../../R'))
print 'RDIR', RDIR
if __name__=='__main__':
    print 'CURRENT_DIR', CURRENT_DIR
    sys.path.append(realpath(join(CURRENT_DIR, '../../')))
    sys.path.append(realpath(join(CURRENT_DIR, '../../../')))
    #sys.path.append(realpath(join(CURRENT_DIR, '../../../../')))
    sys.path.append(realpath(join(CURRENT_DIR, '../../../../R')))
    for p in sys.path: print p
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "tb_website.settings.local")

import unittest
from apps.script_helper.rscript_runner import RScriptRunner
from apps.utils.msg_util import *

class RScriptTest(unittest.TestCase):

    
    def setUp(self):
        os.chdir( RDIR)
        
        self.rscript_file = os.path.join(RDIR, 'test2.r')
        self.test_data_file_ok = os.path.join(RDIR, 'testdata2.csv')
        self.test_data_file_fail = os.path.join(RDIR, 'testdata-fail.csv')

        if not os.path.isfile(self.rscript_file):
            raise Exception('rscript_file not found: %s' % self.rscript_file)

        if not os.path.isfile(self.test_data_file_ok):
            raise Exception('test_data_file not found: %s' % self.test_data_file_ok)

        if not os.path.isfile(self.test_data_file_fail):
            raise Exception('test_data_file not found: %s' % self.test_data_file_fail)
            

    def test1_working_file(self):
        msgt('Test working file: %s' % self.test_data_file_ok)
        runner = RScriptRunner(self.test_data_file_ok)
        runner.run_command_on_file()

        msg(runner.get_formatted_response())

        assert not runner.err_found, True
        assert(runner.get_formatted_response(), {"r":[0.64],"s":[0.36]})


    #@unittest.skip('skipping test2_bad_file')
    def test2_bad_file(self):
        msgt('Test bad file: %s' % self.test_data_file_fail)
        runner = RScriptRunner(self.test_data_file_fail)
        runner.run_command_on_file()
         
        msg(runner.get_err_msgs())
         
        assert runner.err_found, True
         #print (runner.get_formatted_response())
        assert(runner.get_err_msgs(), 'Sorry!  There was was an error running the script.<br />CalledProcess error(1): <br />(file name: testdata-fail.csv)')


    @unittest.skip('skipping test3_nonexistent_file')
    def test3_nonexistent_file(self):
        msgt('Test nonexistent file: blah.csv')
        runner = RScriptRunner('blah.csv')
        runner.run_command_on_file()

        msg(runner.get_err_msgs())

        assert runner.err_found, True
        assert(runner.get_err_msgs(), 'The file was not found: "blah.csv"')


    def test4_non_existent_rscript(self):

        runner = RScriptRunner(self.test_data_file_fail)
        runner.rcommand_base = 'Rscript test2_tis_not_a_real_script.r'
        runner.run_command_on_file()
        #self.rscript_file = os.path.join(RDIR, 'test2.r')
        msg(runner.get_err_msgs())

        assert runner.err_found, True

        assert(runner.get_err_msgs(as_html=True), """Sorry!  There was an error running the script. (5)<br />Error: Fatal error: cannot open<br />file 'test2_tis_not_a_real_script.r': No such file or directory<br /><br />R script path: /Users/rmp553/Documents/iqss-git/PhthisisRavens/R""")


if __name__=='__main__':
    unittest.main()
    #test = RScriptTest()





