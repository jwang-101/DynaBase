"""
Functions to help working with families in dimension 1
over number fields.

AUTHORS:

- Ben Hutz (2023-11): initial version

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


#sagerel -pip install pysha3
import sha3  #adds shake to hashlib
import hashlib  #for shake
from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base


##########################

load("functions/function_dim_1_helpers_generic.py")


###########################################
# Families over Number Fields
###########################################

def normalize_family_NF(F, log_file=sys.stdout):
    """
    put the base field in normalized form.

    TODO: Does this need to return the embedding?
    """
    K=F.base_ring()
    Kbase = F.base_ring().base_ring()
    if Kbase in NumberFields():
        if Kbase != QQ:
            def_poly = Kbase.defining_polynomial()
            def_poly_new, gen = pari(def_poly).polredabs(1)
            def_poly_new = def_poly_new.sage({def_poly.parent().variable_names()[0]:def_poly.parent().gen(0)})
            if def_poly_new != def_poly:
                gen = gen.lift().sage({def_poly.parent().variable_names()[0]:def_poly.parent().gen(0)})
                if emb is None:
                    if Kbase.coerce_embedding() is None:
                        L = NumberField(def_poly_new, 'a', embedding=CC(def_poly_new.roots(ring=CC)[0][0]))
                    else:
                        L = NumberField(def_poly_new, 'a', embedding=gen(Kbase.coerce_embedding()(Kbase.gen())))
                else:
                    L = NumberField(def_poly_new, 'a', embedding=gen(emb(Kbase.gen())))
                to_L = Kbase.hom([gen], L)
                Kbase = Kbase.change_ring(to_L)
                log_file.write('Normalized field:' + str(K) + ' to ' + str(L) + '\n')
            if Kbase.variable_name() != 'a':
                Kbase = Kbase.change_names('a')
                K = K.change_ring(Kbase)
                F = F.change_ring(K)
                log_file.write('Normalized field variable:' + str(K) + '\n')
        N = K.ngens()
        if N > 1 or isinstance(K,MPolynomialRing_base):
            R = PolynomialRing(Kbase,'t',N)
        else:
            R = PolynomialRing(Kbase,'t0')
        if [str(v) for v in K.gens()] != [str(v) for v in R.gens()]:
            F = F.change_ring(R)
        return F
    else: #other field
        raise NotImplementedError('Can only normalize number fields with this function')

def get_sage_family_NF(label, log_file=sys.stdout):
    """
    Given a label and a model name, return the sage dynamical system.

    """
    query={}
    query['label']=label
    my_cursor.execute("""SELECT model_coeffs, num_parameters, base_field_label FROM families_dim_1_NF
        WHERE label=%(label)s
        """, query)
    G = my_cursor.fetchone()
    K = get_sage_field_NF(G['base_field_label'])

    n = G['num_parameters']
    S = PolynomialRing(K,'t',n)

    P = ProjectiveSpace(S,1,'x,y')
    R = P.coordinate_ring()
    x,y = R.gens()
    d = len(G['model_coeffs'][0])-1
    polys = []
    for L in G['model_coeffs']:
        poly = 0
        for i in range(0,d+1):
            poly += x**(d-i)*y**i*S(L[i])
        polys.append(poly)

    return DynamicalSystem(polys, domain=P)


def add_family_NF(F, is_poly=None, num_crit=None, num_aut=None, bool_add_field=False, log_file=sys.stdout, timeout=30):
    """
    Give a family of sage functions F, determine it's label and add it to the database.

    add_field checks to see if the base field is in the database and adds it if it is not


    degree integer,
    num_parameters integer,
    base_field base_field_type,
    sigma_invariants sigma_invariants_type,
    citations integer[],
    model_coeffs varchar[],
    model_resultant varchar,
    is_polynomial boolean,
    num_critical_points integer,
    automorphism_group_cardinality integer
    """
    f = {}
    f['degree'] = int(F.degree())
    f['dimension'] = int(F.codomain().dimension_relative())
    #f['keywords'] = keywords
    F = normalize_family_NF(F)

    base_field = F.base_ring().base_ring()
    f['num_parameters'] = int(base_field.ngens())

    bool, K_id = field_in_database_NF(base_field)
    K_id = K_id
    if not bool:
        if bool_add_field:
            K_id = add_field_NF(base_field, log_file=log_file)
        else:
            log_file.write('Could not add : ' + str(list(F)) + ' because ' + str(base_field) + ' not in database \n')
            raise ValueError("base_field not in database")
    F.normalize_coordinates()

    f['base_field_label'] = K_id
    f['base_field_degree'] = int(base_field.degree())
    #f['base_field_emb'] = int(emb_index)

    f['sigma_invariants.one']=[str(t) for t in F.sigma_invariants(1)]
    sigma_hash = str(hashlib.shake_256(''.join(f['sigma_invariants.one']).encode('utf-8')).hexdigest(digest_length))
    label = str(f['dimension']) +'.'+ str(f['degree']) + '.' + sigma_hash + '.1'

    #see if a conjugate is already in the database
    log_file.write('Searching for functions: ' + label + ': ')
    query = {'degree':f['degree'], 'sigma_invariants.one': f['sigma_invariants.one']}
    my_cursor.execute("""SELECT * FROM families_dim_1_NF
        WHERE degree=%(degree)s AND
            (sigma_invariants).one=%(sigma_invariants.one)s::varchar[]""",query)
    if my_cursor.rowcount != 0:
        F_id = my_cursor.fetchone()['label']
        log_file.write('family already known : ' + str(list(F)) + ' as ' + F_id + '\n')
        return F_id
    # otherwise we'll add the function
    f['label'] = label

    #original model
    f['model_coeffs'] = [get_coefficients(g) for g in F]
    f['model_resultant'] = str(F.resultant())
    log_file.write('model computed: \n')

    my_cursor.execute("""INSERT INTO families_dim_1_NF
        (label, degree, num_parameters, model_coeffs, model_resultant,
         base_field_label, base_field_degree, sigma_invariants.one)
        VALUES
        (%(label)s, %(degree)s, %(num_parameters)s, %(model_coeffs)s, %(model_resultant)s,
         %(base_field_label)s, %(base_field_degree)s, %(sigma_invariants.one)s)
        RETURNING label """,f)
    F_id = my_cursor.fetchone()[0]
    log_file.write('inserted: ' + str(list(F)) + ' as ' + str(F_id) + '\n')

    # update sigmas
    sigma = {}
    sigma['label']=label
    sigma.update({'two':[str(t) for t in F.sigma_invariants(2)]})
    my_cursor.execute("""UPDATE families_dim_1_NF
        SET sigma_invariants.two = %(two)s
        WHERE label=%(label)s
        """,sigma)
    if my_cursor.rowcount == 0:
        log_file.write('function ' + label + 'not found \n')
        raise ValueError("function not found to update")
    else:
        log_file.write('updated ' + str(my_cursor.rowcount) + ' functions for sigma.2\n')

    sigma.update({'three':[str(t) for t in F.sigma_invariants(3)]})
    my_cursor.execute("""UPDATE families_dim_1_NF
        SET sigma_invariants.three = %(three)s
        WHERE label=%(label)s
        """,sigma)
    if my_cursor.rowcount == 0:
        log_file.write('sigma: function ' + label + 'not found \n')
        raise ValueError("sigma: function not found to update")
    else:
        log_file.write('updated ' + str(my_cursor.rowcount) + ' functions for sigma.3\n')

    # update booleans
    my_cursor.execute("""UPDATE families_dim_1_NF
        SET is_polynomial=%s,
            num_critical_points=%s,
            automorphism_group_cardinality=%s
        WHERE label=%s
        """,[is_poly, num_crit, num_aut, label])

    return F_id


def add_citations_family_NF(label, citations, log_file=sys.stdout):
    """
        Add the id of the citations for this family
    """
    #make sure they are ints
    citations = [int(t) for t in citations]

    my_cursor.execute("""SELECT
        citations
         FROM families_dim_1_NF
        WHERE label=%s
        """,[label])
    #merge and sort the new list of citations
    if my_cursor.rowcount == 0:
        new_cites = sorted(citations)
    else:
        cites = my_cursor.fetchone()['citations']
        if cites is None:
            new_cites = sorted(citations)
        else:
            new_cites = sorted(list(set(cites+citations)))
    log_file.write('new citations list for ' + label + ' is ' + str(new_cites) + '\n')
    my_cursor.execute("""UPDATE families_dim_1_NF
        SET citations = %s
        WHERE
            label = %s
        """, [new_cites, label])
