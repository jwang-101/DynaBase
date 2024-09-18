load_attach_path('/mnt/c/Users/glisc/projects/Dynabase')
import sys
sys.path.append("/mnt/c/Users/glisc/projects/Dynabase")

# set up tables
load("connect.py")
load("setup_tables/setup_all.py")

# field data
load("connect.py")
load("sample_data/add_fields.py")

# family data
load("connect.py")
load("sample_data/add_families.py")

# functions data
load("connect.py")
load("sample_data/add_functions_dim_1.py")
