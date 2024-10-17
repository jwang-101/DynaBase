
###################################
###connect to database

load("connect.py")

#returns
# my_session - database connection
# my_cursor - cursor to my_session

path_to_log = "/mnt/c/Users/glisc/projects/DynaBase/families_log.txt"
log_file = open(path_to_log, 'w', 1)

from functions.families_dim_1_helpers_NF import add_family_NF
from functions.families_dim_1_helpers_NF import add_citations_family_NF

#########################


R=PolynomialRing(QQ,1,'c')
c = R.gen()
P=ProjectiveSpace(R,1,'x,y')
x,y = P.gens()
F=DynamicalSystem([x**2+c*y**2,y**2])
family_id = add_family_NF(F, my_cursor, is_poly=True, num_crit=int(2), num_aut=int(1), name='poly_deg_2')

my_cursor.execute("""SELECT
        id
         FROM citations
        WHERE label=%s
        """,['Poonen1998'])
cites = my_cursor.fetchone()
add_citations_family_NF(family_id, cites, my_cursor)

R=PolynomialRing(QQ,1,'c')
c = R.gen()
P=ProjectiveSpace(R,1,'x,y')
x,y = P.gens()
F=DynamicalSystem([x**3+c*y**3,y**3])
family_id = add_family_NF(F, my_cursor, is_poly=True, num_crit=int(2), num_aut=int(1), name='poly_deg_3')


###########################


my_session.commit()
my_session.close()
