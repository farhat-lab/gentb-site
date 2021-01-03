# genTB website

A Django site that facilitates the contribution and analysis of fastQ and VCF files.

## Installation

This website requires python 3.x plus a database with GIS functions.

Create a virtualenv and install deps:

```bash
virtualenv -p python3 pythonenv
./pythonenv/bin/activate
python3 -m pip install --upgrade pip
pip install -r requirements.txt
```

Next, you must create a database if you are using PostgreSQL or MySQL. If you are using the default SQLite, skip this.

```bash
vim ./tb-website/settings/local.py
```

Change settings relating to the database `location`, `username`, and `password`.

Finally, it would be best to make sure GIS libraries are installed (some are only required for SQLite and MySQL). The version of `libgeos` may differ depending on the Ubuntu version.

```bash
sudo apt install binutils libproj-dev gdal-bin libgeos-3.6.2 libsqlite3-mod spatialite
```

Next, you should be able to run the migration:

```bash
python manage.py migrate
```

If successful, everything from here is just a standard Django website.

Start the development server:

```bash
python manage.py runserver
```

The database will be empty, so be sure to populate it with a user account and other information you need.

## Maps Testing

To test the projection, maps and other visualisations you can quickly bring a blank database up with the following commands::

```bash
./manage migrate
./manage load_map_data
./manage loaddata drugs genelocus
./manage loaddata test-genetics test-strains
./manage load_social_data
./manage runserver
```

The populated data will be wrong (not real) but should allow necessary testing.

## MacOS:

Comment out the spatialite_library_path in your local.py file

You will likely need to install the `gdal` library. 

### homebrew

```bash
brew install gdal
```

### macports

```bash
macports install gdal
```

Finally, install the required python library.

```bash
pip install gdal
```
