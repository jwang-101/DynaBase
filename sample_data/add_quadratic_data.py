
from functions.function_dim_1_helpers_NF import model_in_database_NF

from functions.function_dim_1_helpers_NF import add_function_all_NF

from sage.parallel.ncpus import ncpus

from sage.parallel.use_fork import p_iter_fork

###################################
###connect to database

load("connect.py")

#returns
# my_session - database connection
# my_cursor - cursor to my_session

#########################

path_to_log = "/mnt/c/Users/glisc/projects/DynaBase/sample_data/quadratic_log.txt"
log_file = open(path_to_log, 'w', 1)

#x^2 + c
P=ProjectiveSpace(QQ,1,'x,y')
x,y = P.gens()

def parallel_function(F):
    if not model_in_database_NF(F, my_cursor)[0]:
        return add_function_all_NF(F, my_cursor, citations=['Poonen1998'], log_file=log_file)

# Calling possible_periods for each prime in parallel
parallel_data = []
for c in QQ.range_by_height(8): #test bigger numbers here; crash error or timeout error #check if parallelism works without crashing
    F=DynamicalSystem([x**2+c*y**2,y**2]) #polys
    parallel_data.append(((F,), {}))

parallel_iter = p_iter_fork(ncpus()-1, 0) 
parallel_results = list(parallel_iter(parallel_function, parallel_data))

my_session.commit()


log_file.close()

#my_session.close()
