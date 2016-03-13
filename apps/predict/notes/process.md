# 1 - Files uploaded from dropbox links

  - User has specified file type and whether single or pair end

# 2 - Add method to create script, e.g.

## FastQ script

- **analyseNGS.pl**
  - basic script: ```perl analyseNGS.pl (single|pair) (extension) (file source directory)```
  - example, single-ended: ```perl analyseNGS.pl 0 . /gentb/runs/job_00001```
  - example, pair-ended: ```perl analyseNGS.pl 1 . /gentb/runs/job_00001```

## VCF script

- **analyseVCF.pl**
  - basic script: ```perl analyseVCF.pl (file source directory)```
  - example: ```perl analyseVCF.pl /gentb/runs/job_00004```
