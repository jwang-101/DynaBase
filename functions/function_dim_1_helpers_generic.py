"""
Functions to help working with functions in dimension 1
over any field.

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



from sage.rings.fraction_field import is_FractionField
from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general

#sagerel -pip install pysha3
import sha3  #adds shake to hashlib
import hashlib  #for shake


###########################################
# Global Constants
###########################################

#length of SHAKE-256S hash for function label
digest_length = int(4)

###########################################
# General Functions
###########################################

def get_coefficients(F, d=None):
    """
        Get all coefficients of F (including 0) in lexicographic order
    """
    if d is None:
        d = F.degree()
    C = []
    x, y = F.parent().gens()
    for i in range(0,d+1):
        C.append(str(F.coefficient({x:d-i,y:i})))
    return C

def get_post_critical(Fbar):
    """
    Determine the critical and post-critical set of the given map
    """
    set_verbose(None)
    post_crit = set()
    crit = Fbar.critical_points()
    images_needed = copy(crit)
    while len(images_needed) != 0:
        Q = images_needed.pop()
        Q2 = Fbar(Q)
        if Q2 not in post_crit:
           post_crit.add(Q2)
           images_needed.append(Q2)
    return crit, list(post_crit)

def choose_display_model(label, log_file=sys.stdout):
    """
    From the list of computed models set one as display.
    If the model required a field extension over original, then it is not chosen.
    The order of preference is

    chebyshev
    monic centered
    reduced
    original
    """
    #determine which to display
    query={}
    query['label']=label
    my_cursor.execute("""SELECT (original_model).base_field_label FROM functions_dim_1_NF where label = %(label)s""", query)
    original_field_label = my_cursor.fetchone()['base_field_label']
    my_cursor.execute("""SELECT is_chebyshev FROM functions_dim_1_NF where label = %(label)s""",query)
    is_cheby = my_cursor.fetchone()['is_chebyshev']
    if is_cheby:
        query['display_model'] = 'chebyshev'
        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET display_model = %(display_model)s
            WHERE
                label = %(label)s
            """, query)
        return True
    my_cursor.execute("""SELECT is_polynomial FROM functions_dim_1_NF where label = %(label)s""",query)
    is_poly = my_cursor.fetchone()['is_polynomial']
    if is_poly:
        my_cursor.execute("""SELECT (monic_centered).base_field_label FROM functions_dim_1_NF where label = %(label)s""",query)
        mc_field_label = my_cursor.fetchone()['base_field_label']
        if original_field_label == mc_field_label:
            query['display_model'] = 'monic centered'
            my_cursor.execute("""UPDATE functions_dim_1_NF
                SET display_model = %(display_model)s
                WHERE
                    label = %(label)s
                """, query)
            return True
    my_cursor.execute("""SELECT (reduced_model).coeffs FROM functions_dim_1_NF where label = %(label)s""",query)
    my_coeffs = my_cursor.fetchone()['coeffs']
    if not my_coeffs is None:
        my_cursor.execute("""SELECT (reduced_model).base_field_label FROM functions_dim_1_NF where label = %(label)s""",query)
        red_field_label = my_cursor.fetchone()['base_field_label']
        if original_field_label == red_field_label:
            query['display_model'] = 'reduced'
            my_cursor.execute("""UPDATE functions_dim_1_NF
                SET display_model = %(display_model)s
                WHERE
                    label = %(label)s
                """, query)
            return True
    my_cursor.execute("""SELECT (newton_model).coeffs FROM functions_dim_1_NF where label = %(label)s""",query)
    my_coeffs = my_cursor.fetchone()['coeffs']
    if not my_coeffs is None:
        my_cursor.execute("""SELECT (newton_model).base_field_label FROM functions_dim_1_NF where label = %(label)s""",query)
        new_field_label = my_cursor.fetchone()['base_field_label']
        if original_field_label == new_field_label:
            query['display_model'] = 'newton'
            my_cursor.execute("""UPDATE functions_dim_1_NF
                SET display_model = %(display_model)s
                WHERE
                    label = %(label)s
                """, query)
            return True
    query['display_model'] = 'original'
    my_cursor.execute("""UPDATE functions_dim_1_NF
        SET display_model = %(display_model)s
        WHERE
            label = %(label)s
        """, query)
    return True

def graph_to_array(G):
    #graph to array
    n = G.num_verts()
    G.relabel(tuple([t for t in range(n)]))
    E = [-1 for i in range(n)]
    for v in G.edges():
        E[v[0]] = v[1]
    return E

def array_to_graph(E):
    #array to graph
    Ed = []
    for i in range(len(E)):
        Ed.append((i,E[i]))
    return DiGraph(Ed, loops=True)
