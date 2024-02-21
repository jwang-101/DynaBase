# DAD
Database of Arithmetic Dynamics

This database collects arithmetic information on algebraic dynamical systems.
The database in implemented in PSQL, while the computations are done in Python/SageMath.
The code/functionality is written to be run under a SageMath kernel.


# Pre-requisites and Set-up
On an existing 20.02 ubunutu install with release sage built from sage github
(https://github.com/sagemath/sage.git)

# Install PostGreSQL
sudo apt update
sudo apt upgrade
sudo apt install postgresql

# start it
sudo systemctl start postgresql

# check its status
sudo systemctl status postgresql

# From a terminal
# create database and user
sudo su - postgres
psql
postgres=# CREATE USER [username] WITH PASSWORD '[password]';
postgres=# CREATE DATABASE [db name];
postgres=# GRANT ALL PRIVILEGES ON DATABASE [db name] to [user];
postgres=# \q


# Install a GUI client
# Beekeeper community edition
# Get link from readme in: https://github.com/beekeeper-studio/beekeeper-studio
# download .deb and install with
sudo dpkg -i [filename]

# Connecting via python needs the dependencies
sudo apt-get install python3-dev
sudo apt-get install libpq-dev

# Sage needs psql interface
sage -pip install psycopg2

# Sage needs the hash library to create labels
sage -pip install pysha3

# Connect to the DAD github project
# log in to you github account
# navigate to database page
https://github.com/bhutz/DAD.git

Click on "Fork" to create your personal fork of the project
# this creates your own person DAD repository. e.g., https://github.com/alice/DAD.git

# clone the repository
git clone git clone https://github.com/alice/Dynabase.git
cd DAD
git remote -v

# add the remote as the upstream
git remote add upstream https://github.com/bhutz/Dynabase.git
git remote -v

# disable pushes to upstream instead of origin
git remote set-url --push upstream DISABLE

# If you are using SSH for github, then you would do the following instead:
# create local directory for the project
# initialize as a git repository
git init
# add the remotes
git remote add origin git@github.com:alice/Dynabase.git
git remote add upstream git@github.com:bhutz/Dynabase.git
git remote -v
# pull the repository
git pull origin main
# set tracking so you don't have to specify origin/main everytime
git branch --set-upstream-to=origin/main

# to load the .oy files in Sage. It assumes you've added the root of
# the project to the load_attach_path
sage: load_attach_path('[path]')