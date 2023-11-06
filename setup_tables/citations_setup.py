"""
Setup the PostgreSQL table for citations

AUTHORS:

- Ben Hutz (2023-10): initial version

"""

# ****************************************************************************
#       Copyright (C) 2023 Ben Hutz <benjamin.hutz@slu.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# drop table if it already exists
my_cursor.execute("""
    DROP TABLE IF EXISTS citations
""")


############################
# create custom types


#create new types


######################################
# Create Citations table Schema

#should this have things like: base_field_type?
my_cursor.execute("""
CREATE TABLE citations (
    label varchar PRIMARY KEY,
    authors varchar[],
    journal varchar,
    year integer,
    citation varchar,
    mathscinet varchar
  )""")


#'https://mathscinet.ams.org/mathscinet/article?mr=MR2501344'
my_session.commit()

############################################

bibliography=[]
bibliography.append(['Benedetto2009',['Robert Benedetto', 'Ben Dickman', 'Sasha Joseph', 'Ben Krause', 'Dan Rubin', 'Xinwen Zhou'],\
    'Involve', 2009,\
    'Robert Benedetto, Ben Dickman, Sasha Joseph, Ben Krause, Dan Rubin, and Xinwen Zhou. Computing points of small height for cubic polynomials. Involve, 2:37--64, 2009.',\
    'MR2501344'])
bibliography.append(['BFHJY2018', ['Robert Benedetto', 'Xander Faber', 'Benjamin Hutz', 'Jamie Juul', 'Yu Yasafuku'],\
    'Research in Number Theory', 2017,\
    'Robert Benedetto, Xander Faber, Benjamin Hutz, Jamie Juul, Yu Yasafuku. A large arboreal Galois representation for a cubic postcritically finite polynomial.',\
    'MR3736808'])
bibliography.append(['BM2012', ['Nils Bruin', 'Alex Molnar'], 'LMS J. Comput. Math.', 2012,\
    'Nils Bruin, Alex Molnar. Minimal models for rational functions in a dynamical setting. LMS J. Comput. Math. 15 (2012) 400--417.',\
    'MR3015733'])
#bibliography['Canci2010'] = 'Jung Kyu Canci, Rational periodic points for quadratic maps. Annales de l`institut Fourier, Volume 60 (2010) no. 3, p. 953--985.'
#bibliography['Chang2006'] = 'Jianming Chang, Polynomials without repelling periodic point of given period. J. Math. Anal. Appl. 324 (2006) 1--13.'
#bibliography['Dickson1958'] = 'L.E. Dickson. Linear Groups with an Exposition of the Galois Field Theory. New York: Dover, 1958.'
#bibliography['Doyle2014'] = 'John R. Doyle, Xander Faber, and David Krumm. Computation of preperiodic structures for quadratic polynomials over-number fields. New York Journal of Mathematics, 20:507--605, 2014.'
#bibliography['dFH2016'] = 'Joao de Faria, Benjamin Hutz, Automorphism Groups and Invariant Theory on PN. Journal of Algebra and Its Applications. doi.org/10.1142/S0219498818501621.'
#bibliography['FMV2015'] = 'Xander Faber, Michelle Manes, Bianca Viray, Computing conjugating sets and automorphism groups of rational functions. Journal of Algebra, 423:1161--1190.'
#bibliography['GHK2018'] = 'Thomas Gauthier, Benjamin Hutz, Scott Kaschner, Symmetrization of rational maps: arithmetic properties and families of Lattes maps of Pk. Submitted'
#bibliography['Hutz2015'] = 'Benjamin Hutz. Determination of all rational preperiodic points for morphisms of PN. Mathematics of Comptutaion, 84(291):289--308, 2015.'
#bibliography['Hutz2013'] = 'Benjamin Hutz and Patrick Ingram. Numerical evidence for a conjecture of Poonen. Rocky Mountain Journal of Mathematics, 43(1):193--204, 2013.'
#bibliography['Ingram2012'] = 'Patrick Ingram. A finiteness result for post-critically finite polynomials. Int. Math. Res. Not., 3:524--543, 2012.'
#bibliography['Jones2008'] = 'Rafe Jones. The density of prime divisors in the arithmetic dynamics of quadratic polynomials. J. Lond. Math. Soc. (2), 78(2):523--544, 2008.'
#bibliography['JM2011'] = 'Rafe Jones, Michelle Manes. Galois theory of quadratic rational functions. Commentarii Mathematici Helvetici 89(1) 2011.'
#bibliography['Lukas2014'] = 'David Lukas, Michelle Manes, and Diane Yap. A census of quadratic post-critically finite rational functions defined over Q. LMS Journal of Computation and Mathematics, A:314--329, 2014.'
#bibliography['Manes2008'] = 'Michelle Manes. Q-rational cycles for degree-2 rational maps having an automorphism. Proc. London Math. Soc., 96:669--696, 2008.'
#bibliography['Miasnikov2014'] = 'Nikita Miasnikov, Brian Stout, and Phillip Williams. Automorphism loci for the moduli space of rational maps. arxiv:1408.5655, 2014.'
bibliography.append(['Poonen1998', ['Bjorn Poonen'], 'Math. Z.', 1998,\
    'Bjorn Poonen. The classificiation of rational preperiodic points of quadratic polynomials over Q: a refined conjecture. Math. Z., 228(1):11--29, 1998.',\
    'MR1617987'])
#bibliography['Stoll2008'] = 'Michael Stoll. Rational 6-cycles under iteration of quadratic polynomials. London Math. Soc. J. Comput. Math., 11:367--380, 2008.'
#bibliography['WR1994'] = 'Ralph Walde, Paula Russo. Rational Periodic Points of the Quadratic Function Qc(x) = x2 + c. The American Mathematical Monthly, 101(4):318--331 Apr. 1994.'


for row in bibliography:
    my_cursor.execute("""INSERT INTO citations
        (label, authors, journal, year, citation, mathscinet)
        VALUES
        (%s, %s, %s, %s, %s, %s)
        """,row)

my_session.commit()
my_session.close()



