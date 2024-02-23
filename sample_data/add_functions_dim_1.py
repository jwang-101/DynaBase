
from functions.function_dim_1_helpers_NF import model_in_database_NF

from functions.function_dim_1_helpers_NF import add_function_all_NF

###################################
###connect to database

load("connect.py")

#returns
# my_session - database connection
# my_cursor - cursor to my_session

#########################

path_to_log = "/home/ben/DynaBase/sample_data/functions_log.txt"
log_file = open(path_to_log, 'w', 1)

#x^2 + c
P=ProjectiveSpace(QQ,1,'x,y')
x,y = P.gens()
for c in QQ.range_by_height(8):
    F=DynamicalSystem([x**2+c*y**2,y**2]) #polys
    if not model_in_database_NF(F, my_cursor)[0]:
        label = add_function_all_NF(F, my_cursor, citations=['Poonen1998'], log_file=log_file)

my_session.commit()


#all quadratic rational
P=ProjectiveSpace(QQ,1,'x,y')
x,y = P.gens()
B=3
count = 0
for s1 in QQ.range_by_height(B):
    print("s1:",s1, " count:",count)
    my_session.commit()
    for s2 in QQ.range_by_height(B):
        count += 1
        F=DynamicalSystem([2*x**2 + (2-s1)*x*y + (2-s1)*y**2, -x**2+(2+s1)*x*y + (2-s1-s2)*y**2])
        F.normalize_coordinates()
        if F.degree() == 2:
            if not model_in_database_NF(F, my_cursor)[0]:
                label = add_function_all_NF(F, my_cursor, log_file=log_file)

my_session.commit()


#cubic poly
P=ProjectiveSpace(QQ,1,'x,y')
x,y = P.gens()
B=3
count = 0
for s1 in QQ.range_by_height(B):
    print("s1:",s1, " count:",count)
    my_session.commit()
    for s2 in QQ.range_by_height(B):
        # need to fix these (timeout)
        if (s1,s2) not in [(2,0)]:
            count += 1
            F=DynamicalSystem([x**3 + s1*x*y**2 + s2*y**3, y**3])
            if not model_in_database_NF(F, my_cursor)[0]:
                label = add_function_all_NF(F, my_cursor, citations=['Benedetto2009'], log_file=log_file)

my_session.commit()

# twist example
P = ProjectiveSpace(QQ,1,'x,y')
x,y = P.gens()
F = DynamicalSystem([x**2-2*x*y,-2*x*y+y**2])
add_function_all_NF(F, my_cursor, log_file=log_file)
p = 2
q = 4
G = DynamicalSystem([2*p*x*y-3*q*y**2,3*x**2 - p*y**2])
add_function_all_NF(G, my_cursor, log_file=log_file)
p = 6
q = 21
G2 = DynamicalSystem([2*p*x*y-3*q*y**2,3*x**2 - p*y**2])
add_function_all_NF(G2, my_cursor, log_file=log_file)

# newton example
P = ProjectiveSpace(QQ,1,'x,y')
x,y = P.gens()
F = DynamicalSystem([2/3*x**3+x**2*y - 1/3*y**3, x**2*y + 2*x*y**2], domain=P) #Newton
add_function_all_NF(F, my_cursor, log_file=log_file)

#Lattes example
E = EllipticCurve([0,0,0,0,1])
F = P.Lattes_map(E,2)
add_function_all_NF(F, my_cursor, log_file=log_file)


log_file.close()

my_session.commit()
#my_session.close()
