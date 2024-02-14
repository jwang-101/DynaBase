
###################################
###connect to database

load("connect.py")

#returns
# my_session - database connection
# my_cursor - cursor to my_session

path_to_log = "/home/ben/DynaBase/sample_data/families_log.txt"
log_file = open(path_to_log, 'w', 1)

from functions.families_dim_1_helpers_NF import add_family_NF
from functions.families_dim_1_helpers_NF import add_citations_family_NF

#########################


R=PolynomialRing(QQ,1,'c')
c = R.gen()
P=ProjectiveSpace(R,1,'x,y')
x,y = P.gens()
F=DynamicalSystem([x**2+c*y**2,y**2])
label = add_family_NF(F, my_cursor, is_poly=True, num_crit=int(2), num_aut=int(1))

my_cursor.execute("""SELECT
        id
         FROM citations
        WHERE label=%s
        """,['Poonen1998'])
cites = my_cursor.fetchone()
add_citations_family_NF(label, cites, my_cursor)

R=PolynomialRing(QQ,1,'c')
c = R.gen()
P=ProjectiveSpace(R,1,'x,y')
x,y = P.gens()
F=DynamicalSystem([x**3+c*y**3,y**3])
label = add_family_NF(F, my_cursor, is_poly=True, num_crit=int(2), num_aut=int(1))


###########################


my_session.commit()
my_session.close()
