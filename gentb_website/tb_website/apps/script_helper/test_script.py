import sys
import json
import requests
import random

from random import randint


def get_sample_results():

    result_num = randint(1,5)
    if result_num == 1:
        return """{"amk":{"r":[0.9607],"s":[0.0393]},"cap":{"r":[0.3422],"s":[0.6578]},"cip":{"r":[0.625],"s":[0.375]},"emb":{"r":[0.8466],"s":[0.1534]},"gai":{"r":[0.3709],"s":[0.6291]},"inh":{"r":[0.08],"s":[0.92]},"kan":{"r":[0.8174],"s":[0.1826]},"levo":{"r":[0.2274],"s":[0.7726]},"moxi":{"r":[0.0918],"s":[0.9082]},"oflx":{"r":[0.7672],"s":[0.2328]},"pas":{"r":[0.8614],"s":[0.1386]},"pza":{"r":[0.7266],"s":[0.2734]},"rif":{"r":[0.2518],"s":[0.7482]},"str":{"r":[0.9782],"s":[0.0218]}}"""

    if result_num == 2:
        return """{"amk":{"r":[0.7683],"s":[0.2317]},"cap":{"r":[0.5202],"s":[0.4798]},"cip":{"r":[0.6477],"s":[0.3523]},"emb":{"r":[0.09],"s":[0.91]},"gai":{"r":[0.2632],"s":[0.7368]},"inh":{"r":[0.8933],"s":[0.1067]},"kan":{"r":[0.9487],"s":[0.0513]},"levo":{"r":[0.0042],"s":[0.9958]},"moxi":{"r":[0.3219],"s":[0.6781]},"oflx":{"r":[0.0091],"s":[0.9909]},"pas":{"r":[0.5854],"s":[0.4146]},"pza":{"r":[0.163],"s":[0.837]},"rif":{"r":[0.7605],"s":[0.2395]},"str":{"r":[0.8175],"s":[0.1825]}}"""
    if result_num == 3:
        return """{"amk":{"r":[0.2395],"s":[0.7605]},"cap":{"r":[0.6646],"s":[0.3354]},"cip":{"r":[0.887],"s":[0.113]},"emb":{"r":[0.2078],"s":[0.7922]},"gai":{"r":[0.7926],"s":[0.2074]},"inh":{"r":[0.9169],"s":[0.0831]},"kan":{"r":[0.7129],"s":[0.2871]},"levo":{"r":[0.5527],"s":[0.4473]},"moxi":{"r":[0.161],"s":[0.839]},"oflx":{"r":[0.4254],"s":[0.5746]},"pas":{"r":[0.7354],"s":[0.2646]},"pza":{"r":[0.7647],"s":[0.2353]},"rif":{"r":[0.0993],"s":[0.9007]},"str":{"r":[0.527],"s":[0.473]}}"""
    if result_num == 4:
        return """{"amk":{"r":[0.63],"s":[0.37]},"cap":{"r":[0.5182],"s":[0.4818]},"cip":{"r":[0.9586],"s":[0.0414]},"emb":{"r":[0.086],"s":[0.914]},"gai":{"r":[0.5673],"s":[0.4327]},"inh":{"r":[0.4467],"s":[0.5533]},"kan":{"r":[0.7859],"s":[0.2141]},"levo":{"r":[0.9014],"s":[0.0986]},"moxi":{"r":[0.3027],"s":[0.6973]},"oflx":{"r":[0.2766],"s":[0.7234]},"pas":{"r":[0.0793],"s":[0.9207]},"pza":{"r":[0.8069],"s":[0.1931]},"rif":{"r":[0.5958],"s":[0.4042]},"str":{"r":[0.8607],"s":[0.1393]}}"""
    if result_num == 5:
        return """{"amk":{"r":[0.2513],"s":[0.7487]},"cap":{"r":[0.2782],"s":[0.7218]},"cip":{"r":[0.4983],"s":[0.5017]},"emb":{"r":[0.8704],"s":[0.1296]},"gai":{"r":[0.8634],"s":[0.1366]},"inh":{"r":[0.2478],"s":[0.7522]},"kan":{"r":[0.2348],"s":[0.7652]},"levo":{"r":[0.8065],"s":[0.1935]},"moxi":{"r":[0.5807],"s":[0.4193]},"oflx":{"r":[0.5426],"s":[0.4574]},"pas":{"r":[0.0015],"s":[0.9985]},"pza":{"r":[0.8749],"s":[0.1251]},"rif":{"r":[0.409],"s":[0.591]},"str":{"r":[0.7857],"s":[0.2143]}}"""




def send_back_results(val):
    print ('send_back_results')
    print ('---%s---'%val)

    d = json.loads(val)

    json_result_string = '%s' % json.dumps(dict(var='some_val', val=[1,2,4]))

    payload = dict(success=True,
                run_md5=d['run_md5'],   #'afde133c98aac9657e318de2774e687e',
                result_data=get_sample_results()
                )
    url = d['callback_url']  #http://127.0.0.1:8000/predict/my-dataset-run-notification/'
    r = requests.post(url, data=payload)
    print(r.status_code)
    print(r.text)


def send_back_fail_results(val):
    print ('send_back_fail_results')
    print ('---%s---'%val)
    d = json.loads(val)

    json_result_string = '%s' % json.dumps(dict(error_message='The run failed in step 4.  error code: xyz'))
    payload = dict(success=False,
                run_md5=d['run_md5'],   #'afde133c98aac9657e318de2774e687e',
                result_data=json_result_string
                )
    url = d['callback_url']  #http://127.0.0.1:8000/predict/my-dataset-run-notification/'
    r = requests.post(url, data=payload)
    print(r.status_code)
    print(r.text)


if __name__=='__main__':
    cnt = 0
    for val in sys.argv:
        cnt +=1
        if cnt ==2:
            send_back_results(val)
            """
            if random.choice([True, False]) is True:
                print ('true')
                send_back_results(val)
            else:
                print ('false')
                send_back_fail_results(val)
            """