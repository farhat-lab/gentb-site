#!/bin/sh
# This script should be wrapped by another script that
# encloses all of these commands in "scl enable python27"
# and is run by the "plaid" user one time for setup.
# su plaid -l -s /bin/sh -c 'scl enable python27 "path/to/script.sh"'
# See also http://developerblog.redhat.com/2013/02/14/setting-up-django-and-python-2-7-on-red-hat-enterprise-6-the-easy-way/

# When running after initial build, to update libraries, this script can be run as:
# switch to root:
# $ sudo -s
# set python:
# $ scl enable python27 bash
# then in bash as:
# $ sudo ./script3.sh
# The above steps are currently necessary to get Stampy up.  The other libraries can be installed simply with
# $ sudo ./setup3.sh

# This script installs packages required for converting fastq files into
# a tabular, R format.

## Stampy ##
# http://www.well.ox.ac.uk/project-stampy
# Installation notes:
# http://www.well.ox.ac.uk/~gerton/README.txt
# requires python 2.7
# Download URL:
# http://www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz

# To work in Centos7
# echo "Installing Stampy"
# cd /webapps/code/PhthisisRavens/phthisis_website/tb_website/static/Stampy
# tar -xf Stampy-latest.tar
# cd stampy-1.0.27
# make

# To work in Centos6:
# Need to makefile change 
# -       g++ `$(python)-config --ldflags` -pthread -shared $(objs) -o maptools.so
# +       g++ `$(python)-config --ldflags` -L/opt/rh/python27/root/usr/lib64 -pthread -
# 
# >scl enable python27 bash
# >python -V
# >make



## SAMtools ##
# http://biobits.org/samtools_primer.html
# dependencies: GNU curses library, Zlib compression library

# Intentionally using older version of samtools, to avoid indexing bam files failure

echo "Installing SAMtools"
cd /webapps/code/PhthisisRavens/phthisis_website/tb_website/static/SAMtools
tar -xjf samtools-0.1.19.tar.bz2
cd samtools-0.1.19
make install
#export PATH=/webapps/code/PhthisisRavens/phthisis_website/tb_website/static/SAMtools/samtools-1.2:$PATH


## htslib ##
# http://www.htslib.org/download/

echo "Installing htslib"
cd /webapps/code/PhthisisRavens/phthisis_website/tb_website/static/htslib
tar -xjf htslib-1.2.1.tar.bz2
cd htslib-1.2.1
make install
export PATH=/webapps/code/PhthisisRavens/phthisis_website/tb_website/static/htslib/htslib-1.2.1:$PATH


## Platypus ##
# http://www.well.ox.ac.uk/platypus
# Installation notes on same page.
# Download URL:
# http://www.well.ox.ac.uk/bioinformatics/Software/Platypus-latest.tgz
# dependencies: samtools, htslib

echo "Installing Platypus"
cd /webapps/code/PhthisisRavens/phthisis_website/tb_website/static/Platypus
tar -xf Platypus-latest.tar
cd Platypus_0.8.1
./buildPlatypus.sh

# From Project Page:
# "This will take a minute or so, and generate quite a lot of warnings.
# If the build is successful, you will see a message,
#    'Finished building Platypus'.
# Platypus is then ready for variant-calling."

