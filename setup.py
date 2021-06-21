#!/usr/bin/env python3
#
# Copyright (C) 2018 Maha Farhat
#
# Author: Martin Owens
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3.0 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library.
#

import os
import sys

from setuptools import setup
from tb_website.settings import MOD_VERSION, MOD_PACKAGE

# remove MANIFEST. distutils doesn't properly update it when the
# contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

# Move local settings template away
if os.path.isfile('tb_website/settings/local.py'):
    sys.stderr.write("Refusing to do anything while there is a local.py settings file.\n")
    sys.exit(2)

setup(
    name             = MOD_PACKAGE,
    version          = MOD_VERSION,
    description      = 'GenTB DB Module',
    long_description = 'Provide access to the gentb database, not a full website.',
    author           = 'Martin Owens',
    author_email     = 'doctormo@gmail.com',
    url              = 'https://github.com/farhat-lab/gentb-site',
    platforms        = 'linux',
    license          = 'AGPLv3',
    py_modules       = ['gentb',],
    packages         = """
apps
apps.explore
apps.explore.migrations
apps.maps
apps.maps.management
apps.maps.management.commands
apps.maps.migrations
apps.mutations
apps.mutations.management
apps.mutations.management.commands
apps.mutations.migrations
apps.pipeline
apps.pipeline.management
apps.pipeline.management.commands
apps.pipeline.method
apps.pipeline.migrations
apps.predict
apps.predict.management
apps.predict.management.commands
apps.predict.migrations
apps.predict.templatetags
apps.predict.tests
apps.tb_users
apps.tb_users.migrations
apps.uploads
apps.uploads.management
apps.uploads.management.commands
apps.uploads.migrations
apps.uploads.tests
apps.utils
apps.versioner
apps.versioner.templatetags
tb_website
tb_website.formats
tb_website.formats.en
tb_website.management
tb_website.management.commands
tb_website.settings
tb_website.templatetags
    """.split(),
    scripts=[os.path.join('bin', a) for a in os.listdir('bin')],
    install_requires = [
        'logutils==0.3.3',
        'mysqlclient==1.4.4',
        'Django==2.2.24',
        'PyVCF==0.6.8',
        'django-model-utils==3.1.1',
        'requests-file==1.1',
        'requests-ftp==0.3.1',
        'requests>=2.20.0',
        'pytz==2015.4',
    ],
    data_files=[
    ],
    classifiers      = [
      'Operating System :: POSIX :: Linux',
      'Programming Language :: Python',
      'Programming Language :: Python :: 2.7',
    ],
    include_package_data=True,
 )

