"""
https://blogs.dropbox.com/developers/2015/08/new-api-endpoint-shared-link-metadata/
"""

import requests
import json
from dropbox_info import app_key, secret_key, access_token

example_dlink = 'https://www.dropbox.com/s/4tqczonkaeakvua/001.txt?dl=0'
#example_dlink = 'https://www.dropbox.com/sh/vbicdol2e8mn57r/AACeJzBhUpgxTNjj6jHL2UJoa?dl=0'

def show_relevant_files(dropbox_dict):
    assert isinstance(dropbox_dict, dict),\
            "dropbox_json must be a dict, not: {0}".format(type(dropbox_dict))

    if not "is_dir" in dropbox_dict:
        return None

    if dropbox_dict['is_dir'] is False:
        print 'file: {0}'.format(dropbox_dict['path'])


    print dropbox_dict

def run_metadata(dlink):
    global access_token

    params = dict(link=dlink)#  , path='/folder1')
    headers = {'Authorization': 'Bearer {0}'.format(access_token)}
    print headers
    r = requests.post('https://api.dropbox.com/1/metadata/link',
                    data=params,
                    headers=headers
                    )
    print r.status_code
    print r.text

    print json.dumps(json.loads(r.text), indent=4, sort_keys=True)

    show_relevant_files(r.json())
    #rjson = r.json()


if __name__=='__main__':
    run_metadata(example_dlink)

    #link = 'https://www.dropbox.com/sh/r7qb9skx1vc3xxd/AABPxsafyiXtcsJ0Vg20JNloa?dl=0'
    #print """curl -X POST https://api.dropbox.com/1/metadata/link -u {0}:{1} -d #link={2}""".format(app_key, secret_key, link)
