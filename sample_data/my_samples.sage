load_attach_path('/mnt/c/Users/glisc/projects/DynaBase')
import sys
sys.path.append("/mnt/c/Users/glisc/projects/DynaBase")

# set up tables
load("connect.py")
load("setup_tables/setup_all.py")

# quadratic data
load("connect.py")
load("sample_data/add_quadratic_data.py")

