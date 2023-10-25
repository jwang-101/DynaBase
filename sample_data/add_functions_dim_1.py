

###################################
###connect to database

load("connect.py")

#returns
# my_session - database connection
# my_cursor - cursor to my_session

#########################

path_to_log = "/home/ben/Database/sample_data/functions_log.txt"
log_file = open(path_to_log, 'w', 1)

#x^2 + c
for c in QQ.range_by_height(20):
    F=DynamicalSystem([x**2+c*y**2,y**2]) #polys
    label = add_function_all_NF(F)


"""
#all quadratic rational
P=ProjectiveSpace(QQ,1,'x,y')
x,y = P.gens()
B=6
count = 0
for s1 in QQ.range_by_height(B):
    print("s1:",s1, " count:",count)
    my_session.commit()
    for s2 in QQ.range_by_height(B):
        count += 1
        F=DynamicalSystem([2*x**2 + (2-s1)*x*y + (2-s1)*y**2, -x**2+(2+s1)*x*y + (2-s1-s2)*y**2])
        label = add_function_all_NF(F, log_file=log_file)
"""

##cubic poly
#P=ProjectiveSpace(QQ,1,'x,y')
#x,y = P.gens()
#B=5
#count = 0
#for s1 in QQ.range_by_height(B):
#    print("s1:",s1, " count:",count)
#    my_session.commit()
#    for s2 in QQ.range_by_height(B):
#        count += 1
#        F=DynamicalSystem([x**3 + s1*x*y**2 + s2*y**3, y**3])
#        label = add_function_all_NF(F, log_file=log_file)

log_file.close()
#my_session.commit()
#my_session.close()
