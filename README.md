# genTB website

A Django site that facilitates the contribution and analysis of fastQ and VCF files.

## Installation

This website requires python 3.x plus a database with GIS functions.

Create a virtualenv and install deps:

```bash
virtualenv -p python3 pythonenv
. ./pythonenv/bin/activate
python3 -m pip install --upgrade pip
pip install -r requirements.txt
```

Note: PyVCF can now not be installed using newer setuptools, use `pip install setuptools==58` to downgrade setuptools first.

Next, you must create a database if you are using PostgreSQL or MySQL. If you are using the default SQLite, skip this.

```bash
vim ./tb-website/settings/local.py
```

Change settings relating to the database `location`, `username`, and `password`.

Finally, it would be best to make sure GIS libraries are installed (some are only required for SQLite and MySQL). The version of `libgeos` may differ depending on the Ubuntu version.

```bash
sudo apt install binutils libproj-dev gdal-bin libgeos-dev libsqlite3-mod-spatialite
sudo apt install postgresql-14-postgis-3
```

To create a database for postgres, log into psql:

```sql
CREATE DATABASE gentb_local;
CREATE USER gentb WITH PASSWORD 'password';
ALTER ROLE gentb SET client_encoding TO 'utf8';
ALTER ROLE gentb SET default_transaction_isolation TO 'read committed';
ALTER ROLE gentb SET timezone TO 'UTC';
GRANT ALL PRIVILEGES ON DATABASE gentb_local TO gentb;
ALTER DATABASE gentb_local SET timezone TO 'UTC';

\c gentb_local
CREATE EXTENSION IF NOT EXISTS postgis;
```

For running tests against the database:

```sql
CREATE DATABASE test_gentb_local;
GRANT ALL PRIVILEGES ON DATABASE test_gentb_local TO gentb;
ALTER DATABASE test_gentb_local SET timezone TO 'UTC';

\c test_gentb_local
CREATE EXTENSION IF NOT EXISTS postgis;
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
