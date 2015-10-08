#!/bin/sh
echo "Setting up Phthisis Ravens"
# Platform for
# Lightweight
# Applications from
# IQSS
# Data Science
useradd plaid
# EPEL already enabled on HMDC VM
rpm -Uvh http://dl.fedoraproject.org/pub/epel/6Server/x86_64/epel-release-6-8.noarch.rpm
# on HMDC VM, httpd is already installed
#yum install -y python-pip python-devel httpd mod_wsgi ack elinks
echo "Installing Apache"
yum install -y httpd mod_wsgi ack elinks libjpeg-turbo-devel
echo "Installing Python 2.7"
rpm --import http://ftp.scientificlinux.org/linux/scientific/6.4/x86_64/os/RPM-GPG-KEY-sl
yum install -y http://ftp.scientificlinux.org/linux/scientific/6.4/x86_64/external_products/softwarecollections/yum-conf-softwarecollections-1.0-1.el6.noarch.rpm
yum install -y python27
echo "Setting up Django app with Python 2.7"
echo "Installing pip for Python 2.7"
scl enable python27 "easy_install pip"
echo "Install virtualenvwrapper"
scl enable python27 "pip install virtualenvwrapper"

echo "Setup virtualenv directory"
mkdir -p /webapps/virtualenvs
chown plaid /webapps/virtualenvs
mkdir /webapps/code
chown plaid /webapps/code
# in production, use deploy key and clone with git, not https, from GitHub
su plaid -l -s /bin/sh -c 'cd /webapps/code && cp -r /git/PhthisisRavens .'
cp /git/PhthisisRavens/phthisis_website/tb_website/tb_website/settings/template_secret_settings.json /webapps/code/PhthisisRavens/phthisis_website/tb_website/tb_website/settings/template_secret_settings.json
chown plaid:apache /webapps/code/PhthisisRavens/phthisis_website/tb_website/tb_website/settings/template_secret_settings.json
chmod 440 /webapps/code/PhthisisRavens/phthisis_website/tb_website/tb_website/settings/template_secret_settings.json
#
# Create directory for sqlite db
#
echo "Create general data directory"
#
mkdir -p /webapps/data/tb
chown apache /webapps/data/tb
#
echo "Create data directory for sqlite db (apache needs write access"
mkdir -p /webapps/data/tb/sqlite
chown plaid:apache /webapps/data/tb/sqlite
chmod 775 /webapps/data/tb/sqlite
chown plaid:apache /webapps/data/tb/sqlite/tb_website.db3
chmod 660 /webapps/data/tb/sqlite/tb_website.db3

#
# configure apache
#
echo "Configure Apache"
cp /webapps/code/PhthisisRavens/phthisis_website/deploy/vagrant-centos-tb.conf /etc/httpd/conf.d/tb.conf
chown plaid /etc/httpd/conf.d/tb.conf
echo "Create /var/www directory owned by plaid"
mkdir /var/www/tb
chown plaid /var/www/tb
cp /webapps/code/PhthisisRavens/deploy/files/etc/sudoers.d/plaid /etc/sudoers.d
#
# Create directory for uploaded files, writable by apache
#
mkdir -p /webapps/data/tb/tb_uploaded_files
chown apache /webapps/data/tb/tb_uploaded_files
#
# run main setup script as "plaid" user with python 2.7
#
su plaid -l -s /bin/sh -c 'scl enable python27 "/webapps/code/PhthisisRavens/phthisis_website/scripts/setup2.sh"'

#
# run setup of libraries for manipulating fastq files
#
./setup3.sh

#
# fix permissions on database
#
chown plaid:apache /webapps/data/tb/sqlite/tb_website.db3
chmod 660 /webapps/data/tb/sqlite/tb_website.db3

service httpd start
chkconfig httpd on
# on HMDC VM, changed SELinux to "permissive" in /etc/selinux/config
