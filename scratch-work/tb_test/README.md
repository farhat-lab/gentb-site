
### To run:

 - cd into the ```../tb_test``` directory
### Create a virtualenv:

```
module load dev/python/2.7.6
virtualenv venv_test
pip install -r requirements.txt
```

### Run gunicorn:

```
gunicorn --workers=1 -b 0.0.0.0:9001 wsgi:application
```

### check the gentb page:
  - https://gentb.hms.harvard.edu/
