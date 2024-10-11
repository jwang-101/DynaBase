
###################################
###connect to database

load("connect.py")

#returns
# my_session - database connection
# my_cursor - cursor to my_session

path_to_log = "/mnt/c/Users/glisc/projects/DynaBase/fields_log.txt"
log_file = open(path_to_log, 'w', 1)


from fields.field_helpers_FF import add_field_FF


#########################
###########################

# Finite Fields

for p in primes(2,20):
    for n in range(1,5):
        add_field_FF(GF(p**n), my_cursor)


my_session.commit()
my_session.close()
