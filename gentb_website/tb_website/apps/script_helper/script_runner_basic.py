from __future__ import print_function
from subprocess import Popen, PIPE

#
#   Run a script without waiting for the output
#
def run_script(cmd_args, run_test=False):
    print('run_script: ', cmd_args)
    if run_test:
        cmd_args = ['ls', '-la']
    p = Popen(cmd_args)#, stdout=PIPE, stderr=PIPE)

    #print (p.communicate())

"""
python manage.py shell
from apps.script_helper.script_runner_basic import *
p = run_script(None)
p.communicate()


from subprocess import Popen, PIPE
cmd = 'echo {0}'.format('{"file1":"this.tab"}')

cmd = 'python /Users/rmp553/Documents/iqss-git/PhthisisRavens/phthisis_website/tb_website/apps/script_helper/test_script.py {0}'.format('{"file1":"this.tab"}')

p = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
p.communicate()

p = Popen(cmd.split())

{"file1":"this.tab"}
"""