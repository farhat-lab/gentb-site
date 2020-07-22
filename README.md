# genTB website

Basic django site to facilitate the contribution and analysis of fastQ and VCF files.

# Installation

This website requires python 2.7 plus a database with GIS functions.

Create a virtualenv and install deps:

    virtualenv -p python2.7 pythonenv
    ./pythonenv/bin/activate
    pip install -r requirements.txt

Next you must create a database if you are using postgresql or mysql, (if you are using the default sqlite, skip this)

    vim ./tb-website/settings/local.py

Change settings relating to the database location, username and password.

Finally you must make sure GIS libraries are installed (some are only required for sqlite and mysql), version of libgeos may differ depending on Ubuntu version.

    sudo apt install binutils libproj-dev gdal-bin libgeos-3.6.2 libsqlite3-mod-spatialite

Next you should be able to run the migration:

    python manage.py migrate

If this is sucessfull, everything from here is just a normal django website, to run the server:

    python manage.py runserver

The database will be empty, so be sure to populate it with a user account and other information you need.

# Maps Testing

To test the projection, maps and other visualisations you can quickly bring a blank database up with the following commands::

    ./manage migrate
    ./manage load_map_data
    ./manage loaddata drugs genelocus
    ./manage test-genetics test-strains
    ./manage runserver

These data will be wrong (not real), but should allow basic testing.

# macOS:

Comment out the spatialite_library_path in your local.py file

You will likely need to install the gdal library. If using homebrew, follow these steps::

    brew install gdal
    pip3 install --upgrade pip
    pip3 install gdal
