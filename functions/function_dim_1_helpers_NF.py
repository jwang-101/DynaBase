"""
Functions to help working with functions in dimension 1
defined over number fields.

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
# General Functions
###########################################


load("functions/function_dim_1_helpers_generic.py")

##############################################################
#  Functionality for working with functions over number fields
##############################################################

def normalize_function_NF(F, log_file=sys.stdout):
    """
    check that base field is in normalize form. If not, normalize it
    and the dynamical system. Return the normalized dynamical system.

    """
    base_field = F.base_ring()
    if is_FractionField(base_field) or base_field in FunctionFields():
        lcm_den = lcm([t.denominator() for g in F for t in g.coefficients()])
        F.scale_by(lcm_den)
        base_field = base_field.ring()
        F = F.change_ring(base_field)

    if isinstance(base_field, (PolynomialRing_general, MPolynomialRing_base)):
        R = PolynomialRing(base_field.base_ring(), base_field.ngens(), 't')
        phi = base_field.hom(R.gens(), R)
        F = F.change_ring(phi)
        base_field = base_field.base_ring()
        K, phi = normalize_field_NF(base_field, log_file=log_file)
        if not K is base_field:
            #extend to polyring
            phi = F.base_ring().hom(phi, F.base_ring().change_ring(phi.codomain()))
            F = F.change_ring(phi)
            log_file.write('Normalized base field:' + str(F.base_ring()) + ' to ' + str(K) + '\n')
    else:
        K, phi = normalize_field_NF(base_field, log_file=log_file)
        if not K is base_field:
            F = F.change_ring(phi)
            log_file.write('Normalized base field:' + str(base_field) + ' to ' + str(K) + '\n')

    return F, phi


def model_in_database_NF(F, sigma_1=None, conj_fns=None, log_file=sys.stdout):
    """
    Determine if the model F is in the database.

    if it is return 1, label, otherwise return 0, '0'
    """
    if conj_fns is None:
        if sigma_1 is None:
            sigma_1 = F.sigma_invariants(1)
        query = {'degree':int(F.degree()),\
            'sigma_invariants.one': [str(t) for t in sigma_1]}
        my_cursor.execute("""SELECT
            label,
            (original_model).coeffs,(original_model).base_field_label,
            (monic_centered).coeffs,(monic_centered).base_field_label,
            (chebyshev_model).coeffs,(chebyshev_model).base_field_label,
            (reduced_model).coeffs,(reduced_model).base_field_label,
            (newton_model).coeffs,(newton_model).base_field_label
             FROM functions_dim_1_NF
            WHERE degree=%(degree)s AND
                (sigma_invariants).one=%(sigma_invariants.one)s::varchar[]""",query)
        conj_fns = my_cursor.fetchall()

    model_names = ['original_model', 'monic_centered', 'chebyshev_model', 'reduced_model', 'newton_model']
    #field_label, emb_index, field_id = get_field_label(F.base_ring(), log_file=log_file)
    F_coeffs = [get_coefficients(g) for g in F]
    bool, K_id = field_in_database_NF(F.base_ring(), log_file=log_file)
    if not bool:
        return 0, '0'
    for g in conj_fns:
        for i in range(5):
            if g[2*i+1] == F_coeffs and g[2*i+2] == str(K_id):
                return 1, g[0]
    return 0, '0'


def get_sage_func_NF(label, model_name, log_file=sys.stdout):
    """
    Given a label and a model name, return the sage dynamical system.

    """
    query={}
    query['label']=label
    if model_name == 'original':
        my_cursor.execute("""SELECT (original_model).coeffs, (original_model).base_field_label FROM functions_dim_1_NF
            WHERE label=%(label)s
            """, query)
    elif model_name == 'monic_centered':
        my_cursor.execute("""SELECT (monic_centered).coeffs, (monic_centered).base_field_label FROM functions_dim_1_NF
            WHERE label=%(label)s
            """, query)
    elif model_name == 'chebyshev':
        my_cursor.execute("""SELECT (chebyshev_model).coeffs, (chebyshev_model).base_field_label FROM functions_dim_1_NF
            WHERE label=%(label)s
            """, query)
    elif model_name == 'reduced':
        my_cursor.execute("""SELECT (reduced_model).coeffs, (reduced_model).base_field_label FROM functions_dim_1_NF
            WHERE label=%(label)s
            """, query)
    elif model_name == 'newton':
        my_cursor.execute("""SELECT (newton_model).coeffs, (newton_model).base_field_label FROM functions_dim_1_NF
            WHERE label=%(label)s
            """, query)
    G = my_cursor.fetchone()

    K = get_sage_field_NF(G['base_field_label'])
    P = ProjectiveSpace(K,1,'x,y')
    R = P.coordinate_ring()
    x,y = R.gens()
    d = len(G['coeffs'][0])-1
    polys = []
    for L in G['coeffs']:
        poly = 0
        for i in range(0,d+1):
            poly += x**(d-i)*y**i*Integer(L[i])
        polys.append(poly)

    return DynamicalSystem(polys, domain=P)

def check_conjugates_NF(F,G, normalize_base=False, log_file=sys.stdout):
    """
    F,G are two sage models with same degree, dimension

    returns 0 for different conj class
    1 for rationally conjugate
    2 for rational twists

    """
    if normalize_base:
        F, phiF = normalize_function_NF(F)
        G, phiG = normalize_function_NF(G)
    Kf = F.base_ring()
    Kg = G.base_ring()
    try:
        Fbar = F.change_ring(QQbar)
    except ValueError:
        Fbar = F.change_ring(F.base_ring().embeddings(QQbar)[0])
    try:
        Gbar = G.change_ring(QQbar)
    except ValueError:
        Gbar = G.change_ring(G.base_ring().embeddings(QQbar)[0])
    CS = Fbar.conjugating_set(Gbar)

    if len(CS) == 0:
        return 0
    if Kf != Kg:
        return 1 #either conjugate or not, but not rational twists since defined over different fields
    #so we assume they have the same base field so we are looking for twists
    for m in CS:
        K, E, phi = number_field_elements_from_algebraics([t for r in m for t in r])
        K, psi = normalize_field_NF(K, log_file=log_file)
        if K.is_subring(Kf):
            return 1
    return 2



def conj_in_database_NF(F, conj_fns=None, log_file=sys.stdout, timeout=30):
    """
    Determine if F is conjugate to a model already in the database.
    This includes the identity conjugation.

    This can either be conjugate over the base field or a twist
    (conjugate over an extension of the base field).

    Need to be careful that max_sigma is not larger than what is computed in the database.
    If it is, it will always say the function is new.

    Returns:
    0,[]  if not there.
    1,[g] if rationally conjugate to g
    2,[...] if a twist of ...
    """
    base_field = F.base_ring()
    query = {}
    query['degree'] = int(F.degree())
    sigma_1 = [str(t) for t in F.sigma_invariants(1)]
    query['sigma_invariants.one'] = sigma_1
    if conj_fns is None:
        my_cursor.execute("""SELECT * FROM functions_dim_1_NF
        WHERE degree=%(degree)s AND
            (sigma_invariants).one=%(sigma_invariants.one)s::varchar[]""",query)
        conj_fns = my_cursor.fetchall()
    if len(conj_fns) == 0:
        #nothing with the same sigma_1
        return 0, []

    bool, g_label = model_in_database_NF(F, conj_fns=conj_fns, log_file=log_file)
    if bool:
        return 1, g_label

    #now check for twists/conjugates
    #sigma_2
    sigma_2 = [str(t) for t in F.sigma_invariants(2)]
    query['sigma_invariants.two'] = sigma_2
    my_cursor.execute("""SELECT * FROM functions_dim_1_NF
        WHERE degree=%(degree)s AND
            (sigma_invariants).one=%(sigma_invariants.one)s::varchar[] AND
            ((sigma_invariants).two=%(sigma_invariants.two)s::varchar[] OR
             (sigma_invariants).two IS NULL)""",query)
    conj_fns = my_cursor.fetchall()
    if len(conj_fns) == 0:
        #nothing with the same sigma_1,sigma_2
        return 0, []
    # TODO should we go to sigma_3??

    try:
        if timeout != 0:
            alarm(timeout)
        for g in conj_fns:
            # TODO allow other models?
            twist_val = check_conjugates_NF(get_sage_func_NF(g[0], 'original'),F)
            if twist_val == 1:
                log_file.write('already there: conjugate found in db: ' + str(list(F)) + ' as ' + g[0] + '\n')
                cancel_alarm()
                return 1, [g['label']]
            elif twist_val == 2:
                if g['rational_twists'] is None:
                    twist_list = [g['label']]
                else:
                    twist_list = g['rational_twists'] + [g['label']]
                    twist_list.sort()
                log_file.write('has twist in db: ' + str(list(F)) + ' as ' + str(twist_list) + '\n')
                cancel_alarm()
                return 2, twist_list
        #not in database
        cancel_alarm()
        return 0, []

    except AlarmInterrupt:
        log_file.write('timeout: func_in_db: ' + str(timeout) + ':' + polys['val'] + '\n')
        raise


def add_function_NF(F, bool_add_field=False, log_file=sys.stdout, timeout=30):
    """
    Give a sage function F, determine it's label and add it to the database.

    add_field checks to see if the field is in the database and adds it if it is not


    'degree'
    'label'
    'base_field_label'
    'base_field_degree'
    'twists' - and updates twists field for the twists of F
    'models.original'
        'coeffs'
        'resultant'
        'bad_primes'
        'height'
        'base_field_label'
        'conjugation_from_original'
            'val'
            'base_field_label'
    'display_model' = 'original'
    'sigma_invariants.one'
    'citations'

    """
    f = {}
    f['degree'] = int(F.degree())
    f['dimension'] = int(F.codomain().dimension_relative())
    #f['keywords'] = keywords

    base_field = F.base_ring()
    #if isinstance(base_field, (PolynomialRing_general, MPolynomialRing_base))\
    #  or base_field in FunctionFields() or is_FractionField(base_field):
    #    f['num_parameters'] = int(base_field.ngens())
    #    F, phi = normalize_function(F, log_file=log_file)
    #    base_field = F.base_ring().base_ring()
    #else:
    #    f['num_parameters'] = int(0)
    #    F, phi = normalize_function(F, log_file=log_file)
    #    base_field = F.base_ring()

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
    label = str(f['dimension']) +'.'+ str(f['degree']) + '.' + sigma_hash + '.'

    #see if a conjugate is already in the database
    log_file.write('Searching for functions: ' + label + ': ')
    query = {'degree':f['degree'], 'sigma_invariants.one': f['sigma_invariants.one']}
    my_cursor.execute("""SELECT * FROM functions_dim_1_NF
        WHERE degree=%(degree)s AND
            (sigma_invariants).one=%(sigma_invariants.one)s::varchar[]""",query)
    conj_fns = my_cursor.fetchall()
    m = len(conj_fns)
    found, L = conj_in_database_NF(F, conj_fns=conj_fns, log_file=log_file, timeout=timeout)

    if found == 1:
        #function or rational conjugate already there
        print(L)
        log_file.write('function already known for : ' + L[0] + '\n')
        return L[0]
    # otherwise we'll add the function
    # assume none are conjugate if we get to here
    f['label'] = label + str(m+1)

    #original model
    f['original_model.coeffs'] = [get_coefficients(g) for g in F]
    f['original_model.resultant'] = str(F.resultant())
    if F.base_ring().degree() == 1:
        bad_primes = F.primes_of_bad_reduction()
    else:
        bad_primes = list(set([p.norm() for p in F.primes_of_bad_reduction()])) #remove duplicates
        bad_primes.sort()
    f['original_model.bad_primes'] = [int(p) for p in bad_primes]
    f['original_model.height'] = float(F.global_height())
    f['original_model.base_field_label'] = f['base_field_label']
    #f['original_model.base_field.degree'] = f['base_field.degree']
    #models['original'].update({'base_field_emb': int(emb_index)})
    M = MatrixSpace(F.base_ring(), F.domain().dimension_relative()+1, F.domain().dimension_relative()+1).one()
    #conjugation to original model
    f['original_model.conjugation_from_original'] = [str(t) for r in M for t in r]
    f['original_model.conjugation_from_original_base_field_label'] = f['base_field_label']

    f['display_model'] = 'original'
    log_file.write('original computed: \n')

    my_cursor.execute("""INSERT INTO functions_dim_1_NF
        (label, degree, base_field_label, base_field_degree,
         sigma_invariants.one, citations, family,
         original_model.coeffs, original_model.resultant, original_model.bad_primes,
         original_model.height, original_model.base_field_label,
         original_model.conjugation_from_original,
         original_model.conjugation_from_original_base_field_label, display_model)
        VALUES
        (%(label)s, %(degree)s, %(base_field_label)s, %(base_field_degree)s,
         %(sigma_invariants.one)s, %(citations)s, %(family)s,
         %(original_model.coeffs)s, %(original_model.resultant)s, %(original_model.bad_primes)s,
         %(original_model.height)s, %(original_model.base_field_label)s,
         %(original_model.conjugation_from_original)s,
         %(original_model.conjugation_from_original_base_field_label)s, %(display_model)s)
        RETURNING label """,f)
    F_id = my_cursor.fetchone()[0]
    log_file.write('inserted: ' + str(list(F)) + ' as ' + str(F_id) + '\n')

    if found == 2:
        #There is one or more twists in database
        f['twists'] = L
        log_file.write('twist list: ' + str(L) + ' for ' + f['label'] + '\n')
        my_cursor.execute("""UPDATE functions_dim_1_NF
                    SET rational_twists = %(twists)s
                    WHERE label=%(label)s
                    """,f)
        for twist_label in L:
            my_cursor.execute("""SELECT rational_twists FROM functions_dim_1_NF
                WHERE label=%s """,[twist_label])
            h = my_cursor.fetchone()
            twist_list = copy(L)
            twist_list.remove(twist_label)
            twist_list.append(f['label'])
            twist_list.sort()
            my_cursor.execute("""UPDATE functions_dim_1_NF
                    SET rational_twists = %s
                    WHERE label=%s
                    """,[twist_list, twist_label])
            log_file.write('updating twist list for: ' + twist_label + '\n')
    return F_id


def add_sigma_inv_NF(label, model_name='original', start=2, end=3, log_file=sys.stdout, timeout=30):
    """
    Compute the sigma invariants for the given function specified either by label.

    model specifies which model name to use for computation

    action can be add/replace

    #TODO add check to see if the sigma invariants are already known for 'add'
    """
    if timeout != 0:
        alarm(timeout)
    log_file.write('computing sigmas for : ' + label + ' from ' + str(start) + ' to ' + str(end) + '\n')
    try:
        F = get_sage_func_NF(label,model_name, log_file=log_file)
        sigma = {'label':label}
        for k in range(start, end+1):
            if k == 1:
                raise ValueError("Cannot change sigma.one")
            if k == 2:
                sigma.update({'two':[str(t) for t in F.sigma_invariants(k)]})
                my_cursor.execute("""UPDATE functions_dim_1_NF
                    SET sigma_invariants.two = %(two)s
                    WHERE label=%(label)s
                    """,sigma)
                if my_cursor.rowcount == 0:
                    log_file.write('function ' + label + 'not found \n')
                    raise ValueError("function not found to update")
                else:
                    log_file.write('updated ' + str(my_cursor.rowcount) + ' functions for sigma.' + str(k) + '\n')
            elif k == 3:
                sigma.update({'three':[str(t) for t in F.sigma_invariants(k)]})
                my_cursor.execute("""UPDATE functions_dim_1_NF
                    SET sigma_invariants.three = %(three)s
                    WHERE label=%(label)s
                    """,sigma)
                if my_cursor.rowcount == 0:
                    log_file.write('sigma: function ' + label + 'not found \n')
                    raise ValueError("sigma: function not found to update")
                else:
                    log_file.write('updated ' + str(my_cursor.rowcount) + ' functions for sigma.' + str(k) + '\n')
            else:
                raise NotImplementedError("only upto k=3")

        cancel_alarm()
        return True
    except (AlarmInterrupt, KeyboardInterrupt):
        log_file.write('timeout: sigma inv for:' + str(timeout) + ':' + label + '\n')
    except:
        log_file.write('failure: sigma inv for:' + label + '\n')
        raise

    cancel_alarm()
    return False

def add_is_pcf(label=None, model_name='original', log_file=sys.stdout, timeout=30):
    """
    Determine if the given function (identified by label) is postcritically finite.

    'is_pcf'

    """
    if timeout != 0:
        alarm(timeout)
    try:
        log_file.write('computing is_pcf for : ' + label + '\n')
        F = get_sage_func_NF(label, model_name, log_file=log_file)
        pcf={'label':label}
        try:
            is_pcf = F.is_postcritically_finite()
        except ValueError:
            is_pcf = F.is_postcritically_finite(embedding=F.base_ring().embeddings(QQbar)[0])
        pcf['is_pcf']=is_pcf
        my_cursor.execute("""UPDATE functions_dim_1_NF
                    SET is_pcf = %(is_pcf)s
                    WHERE label=%(label)s
                    """,pcf)

        if my_cursor.rowcount == 0:
            log_file.write('PCF: function ' + label + 'not found \n')
            raise ValueError("PCF: function not found to update")
        else:
            log_file.write('updated ' + str(my_cursor.rowcount) + ' functions for is_pcf \n')
        log_file.write('is_pcf finished\n')
        cancel_alarm()
        return True
    except AlarmInterrupt:
        log_file.write('timeout: is_pcf for:' + str(timeout) + ':' + label + '\n')
    except:
        log_file.write('failure: is_pcf for:' + label + '\n')
        raise

    cancel_alarm()
    return False


def add_critical_portrait(label, model_name='original', log_file=sys.stdout, timeout=30):
    """
    If the function is pcf create the critical point portrait.

    computes:
      critical_portrait_cardinality
      post_critical_cardinality
      critical_portrait_components
      critical_portrait_structure
      critical_portrait_graph_id

    action can be add/replace
    """
    if timeout != 0:
        alarm(timeout)

    log_file.write('computing critical portrait for:' + label + '\n')
    query={}
    query['label']=label
    my_cursor.execute("""SELECT is_pcf FROM functions_dim_1_NF
            WHERE label=%(label)s""",query)
    is_pcf = my_cursor.fetchone()['is_pcf']
    if not is_pcf:
        log_file.write('critical portrait: not pcf:' + label + '\n')
        return True
        #raise ValueError('Function ' + label + ' not pcf')

    cpp = []
    try:
        F = get_sage_func_NF(label, model_name=model_name, log_file=log_file)
        g = F.critical_point_portrait()
        #TODO add graph_id information!!!
        ##query['critical_portrait_graph_id'] = TODO
        try:
            Fbar = F.change_ring(QQbar)
        except ValueError:
            Fbar = F.change_ring(F.base_ring().embeddings(QQbar)[0])
        for q in Fbar.critical_points():
            cpp.append(q.is_preperiodic(Fbar,return_period=True))
        query['critical_portrait_cardinality'] = int(len(g.vertices()))
        query['post_critical_cardinality'] = len([t for t in g.in_degree() if t != 0])
        query['critical_portrait_components'] = [int(t) for t in g.connected_components_sizes()]
        query['critical_portrait_structure'] = cpp

    except PariError:
        log_file.write('pari error create critical portrait:' +  label + '\n')
        pass #ran out of memory pari.allocatemem(##) for more
        #or normal_form was not cooercible?
    except AlarmInterrupt:
        log_file.write('timeout: create critical portrait:' + str(timeout) + ':' + label + '\n')

    if cpp != []:
        ##TODO add graph_id
        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET critical_portrait_cardinality = %(critical_portrait_cardinality)s,
                post_critical_cardinality = %(post_critical_cardinality)s,
                critical_portrait_components = %(critical_portrait_components)s,
                critical_portrait_structure = %(critical_portrait_components)s
            WHERE label=%(label)s
            """,query)

        cancel_alarm()
        log_file.write('critical portrait added:' + label + '\n')
        return True
    return False


def add_automorphism_group_NF(label, model_name='original', log_file=sys.stdout, timeout=30):
    """
    Find the automorphisms group.

    """
    query={}
    query['label']=label
    if timeout != 0:
        alarm(timeout)
    log_file.write('starting aut group for:' + label + '\n')
    try:
        F = get_sage_func_NF(label,model_name=model_name, log_file=log_file)
        try:
            Fbar = F.change_ring(QQbar)
        except ValueError:
            Fbar = F.change_ring(F.base_ring().embeddings(QQbar)[0])

        aut = Fbar.automorphism_group()
        query['automorphism_group_cardinality'] = int(len(aut))
        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET automorphism_group_cardinality = %(automorphism_group_cardinality)s
            WHERE label=%(label)s
            """,query)
        #TODO check rowcount for success
        cancel_alarm()
        log_file.write('aut group computed for:' + label + '\n')
        return True
    except AlarmInterrupt:
        log_file.write('timeout: aut group for:' + str(timeout) + ':' + f['label'] + '\n')
    except:
        log_file.write('failure: aut group for:' + label + '\n')
        raise
    cancel_alarm()
    return False

def identify_graph(G,f,log_file=sys.stdout):
    """
    determine if the digraph is already in the table and returns it's graph_id

    If it is not in the table, then add it.

    G is the graph of perperiodic points
    f is the fuction
    """
    graph_data = {}
    graph_data['cardinality'] = len(G.vertices())
    graph_data['preperiodic_components'] = G.connected_components_sizes()
    graph_data['num_components'] = len(G.connected_components())
    periodic = set()
    for T in G.all_simple_cycles():
        periodic=periodic.union(set(T))
    periodic_cycles = []
    for T in G.connected_components():
        num_periodic=0
        for t in T:
            if t in periodic:
                num_periodic+=1
        periodic_cycles.append(int(num_periodic))
    graph_data['periodic_cycles'] = periodic_cycles
    max_tail = 0
    t = 0
    for Q in G.vertices():
        if G.in_degree(Q) == 0:
            t = 0
            while Q not in periodic:
                t += 1
                Q = f(Q)
            if t > max_tail:
                max_tail = t
    graph_data['max_tail'] = int(max_tail)
    #check whether the graph is in the table
    my_cursor.execute("""SELECT
         graph_id, edges
         FROM graphs_dim_1_NF
        WHERE cardinality=%(cardinality)s AND
            periodic_cycles = %(periodic_cycles)s AND
            num_components = %(num_components)s AND
            preperiodic_components = %(preperiodic_components)s AND
            max_tail = %(max_tail)s
        """,graph_data)
    # Check for isomorphic graphs
    for row in my_cursor.fetchall():
        G_graph = array_to_graph(row['edges'])
        if G_graph.is_isomorphic(G):
            log_file.write('preperiodic graph already in table: ' + str(row['graph_id']) + '\n')
            return row['graph_id']
    # the graph is not in the table, so add it
    # edges relabels to graph verticies so this has to be done last
    graph_data['edges'] = graph_to_array(G)
    my_cursor.execute("""INSERT INTO graphs_dim_1_NF
        (cardinality, edges, num_components, periodic_cycles,
         preperiodic_components, max_tail)
        VALUES
        (%(cardinality)s, %(edges)s, %(num_components)s, %(periodic_cycles)s,
         %(preperiodic_components)s, %(max_tail)s)
        RETURNING graph_id """,graph_data)
    log_file.write('adding preperiodic graph to table: ' + str(graph_data['edges']) + '\n')
    return my_cursor.fetchone()[0]

def add_rational_preperiodic_points_NF(label, model_name='original', field_label=None, log_file=sys.stdout, timeout=30):
    """
    Find the rational preperiodic points and add

    into functions_dim_1_NF table:
        rational_preriodic_graph_id varchar,

    into rational_preperiodic_dim_1_NF table:
        rational_periodic_points varchar[][2],

    into graphs_dim_1_NF table: #if not already there
        cardinality
        edges
        num_components
        periodic_cycles
        preperiodic components
        max_tail
    """
    if timeout != 0:
        alarm(timeout)

    log_file.write('starting rational preperiodic points for:' + label + '\n')
    try:
        graph_data = {}
        preperiodic_data = {}
        query = {}
        query['label']=label
        F = get_sage_func_NF(label, model_name=model_name, log_file=log_file)
        if field_label is None:
            K = F.base_ring()
        else:
            K = get_sage_field_NF(field_label)
            # TODO: this may have embedding issues in some cases
            FK.change_ring(K)
        if timeout != 0:
            alarm(timeout)
        if K.degree() > 4:
            preper = F.rational_preperiodic_graph(prime_bound=[1,5], lifting_prime=5)
        elif K.degree() > 3:
            preper = F.rational_preperiodic_graph(prime_bound=[1,15], lifting_prime=7)
        else:
            preper = F.rational_preperiodic_graph()
        ## TODO: add graph_id

        preperiodic_data['function_label'] = label
        preperiodic_data['base_field_label'] = field_label
        periodic = []
        for c in preper.all_simple_cycles():
            periodic.append([str(t) for t in c[0]])
        preperiodic_data['rational_periodic_points'] = periodic
        # one point per connected component
        # TODO: needs to be the same order as the components are listed in the graph table

        #identify graph and add if necessary
        graph_id = identify_graph(preper, F, log_file=log_file)
        preperiodic_data['graph_id'] = graph_id

        #TODO check that it isn't already there
        my_cursor.execute("""INSERT INTO rational_preperiodic_dim_1_NF
            (function_label, base_field_label, rational_periodic_points, graph_id)
            VALUES
            (%(function_label)s, %(base_field_label)s, %(rational_periodic_points)s, %(graph_id)s)
            RETURNING id """,preperiodic_data)

        #TODO check rowcount for success
        log_file.write('rational preperiodic points computed for:' + label + '\n')
        cancel_alarm()
        return True
    except AlarmInterrupt:
        log_file.write('timeout: preperiodic points for:' + str(timeout) + ':' + label + '\n')
    except:
        log_file.write('failure: preperiodic points for:' + label + '\n')
        raise

    cancel_alarm()
    return False

def add_reduced_model_NF(label, model_name='original', log_file=sys.stdout, timeout=30):
    """
    Compute the reduced model

    coeffs      varchar[],
    resultant   varchar,
    bad_primes  integer[],
    height      double precision,
    base_field_label varchar,
    conjugation_from_original varchar[],
    conjugation_from_original_base_field_label varchar(%s)

    """
    if timeout != 0:
        alarm(timeout)

    query={}
    query['label']=label
    log_file.write('Computing reduced model for:' + label + '\n')
    try:
        F = get_sage_func_NF(label, model_name=model_name, log_file=log_file)
        try:
            log_file.write('trying reduced with dynatomic \n')
            g, M = F.reduced_form(smallest_coeffs=True, dynatomic=True,\
                                  return_conjugation=True, start_n=1)
        except:
            log_file.write('trying reduced with periodic \n')
            g, M = F.reduced_form(smallest_coeffs=True, dynatomic=False,\
                                   return_conjugation=True, start_n=1)
        if M.base_ring() == ZZ:
            M = M.change_ring(QQ)
        log_file.write('reduced form computed \n')
        if g.base_ring().degree() == 1:
            bad_primes = g.primes_of_bad_reduction()
        else:
            bad_primes = list(set([p.norm() for p in g.primes_of_bad_reduction()])) #remove duplicates
            bad_primes.sort()

        g.normalize_coordinates()
        #original model
        bool, K_id = field_in_database_NF(F.base_ring())
        assert(bool)
        query['reduced_model.coeffs'] = [get_coefficients(g) for g in F]
        query['reduced_model.resultant'] = str(F.resultant())
        query['reduced_model.bad_primes'] = [int(p) for p in bad_primes]
        query['reduced_model.height'] = float(F.global_height())
        query['reduced_model.base_field_label'] = K_id
        #query['reduced_model.base_field.degree'] = int(F.base_ring().degree())
        #models['original'].update({'base_field_emb': int(emb_index)})
        #conjugation to original model
        query['reduced_model.conjugation_from_original'] = [str(t) for r in M for t in r]
        #can these fields be different?
        assert(M.base_ring()==F.base_ring())
        query['reduced_model.conjugation_from_original_base_field_label'] = K_id

        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET reduced_model.coeffs = %(reduced_model.coeffs)s,
                reduced_model.resultant = %(reduced_model.resultant)s,
                reduced_model.bad_primes = %(reduced_model.bad_primes)s,
                reduced_model.height = %(reduced_model.height)s,
                reduced_model.base_field_label = %(reduced_model.base_field_label)s,
                reduced_model.conjugation_from_original = %(reduced_model.conjugation_from_original)s,
                reduced_model.conjugation_from_original_base_field_label = %(reduced_model.conjugation_from_original_base_field_label)s
            WHERE
                label = %(label)s
            """, query)
        log_file.write('reduced model computed: ' + label + '\n')
        cancel_alarm()
        #TODO check rowcount
        return True

    except AlarmInterrupt:
        log_file.write('reduced model timeout: ' + str(timeout) + ':' + label + '\n')
    except:
        log_file.write('reduced model failure: ' + label + '\n')
        raise
    cancel_alarm()
    return False

def add_is_polynomial_NF(label, model_name='original', log_file=sys.stdout, timeout=30):
    """
    Determine if the map is a polynomial map (totally ramified fixed point)

    'is_polynomial'
    """
    if timeout != 0:
        alarm(timeout)
    log_file.write('starting is_polynomial for:' + label + '\n')

    try:
        query={}
        query['label']=label
        F = get_sage_func_NF(label, model_name=model_name, log_file=log_file)
        is_poly = F.is_polynomial()

        query['is_polynomial'] = is_poly
        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET is_polynomial = %(is_polynomial)s
            WHERE label=%(label)s
            """,query)
        #TODO check rowcount for success
        log_file.write('is polynomial computed: ' + label + '\n')
        cancel_alarm()
        return True

    except AlarmInterrupt:
        log_file.write('is_polynomial timeout: ' + str(timeout) + ':' + label + '\n')
    except:
        log_file.write('is_polynomial failure: ' + label + '\n')
        raise
    cancel_alarm()
    return False

def add_monic_centered_model_NF(label, model_name='original', log_file=sys.stdout, timeout=30):
    """
    Compute the monic centered model.

    If it hasn't already been checked to be, add_is_polynomial_NF is called.

    coeffs      varchar[],
    resultant   varchar,
    bad_primes  integer[],
    height      double precision,
    base_field_label varchar,
    conjugation_from_original varchar[],
    conjugation_from_original_base_field_label varchar(%s)

    #TODO Note that this has to start from original or the conjugation is wrong

    """
    if timeout != 0:
        alarm(timeout)
    query={}
    query['label'] = label

    my_cursor.execute("""SELECT is_polynomial FROM functions_dim_1_NF where label = %(label)s""",query)
    is_poly= my_cursor.fetchone()['is_polynomial']
    if is_poly is None:
        add_is_polynomial_NF(label, model_name=model_name, log_file=log_file, timeout=timeout)
        my_cursor.execute("""SELECT is_polynomial FROM functions_dim_1_NF where label = %(label)s""",query)
        is_poly= my_cursor.fetchone()['is_polynomial']
    if not is_poly:
        cancel_alarm()
        #no monic centered model if not a polynomial
        return True

    monic_centered = {}
    try:
        F = get_sage_func_NF(label, model_name=model_name, log_file=log_file)
        N = F.domain().dimension()
        #todo fix normal form so this does not have to be special cased
        G,M,phi = F.normal_form(return_conjugation=True)
        #base field may have changed so normalize again
        G, phi = normalize_function_NF(G, log_file=log_file)
        L = G.base_ring()
        M = matrix(L, N+1, N+1, [phi(t) for r in M for t in r])

        G.normalize_coordinates()
        res = G.resultant()
        #correction to normal form
        G.scale_by(1/G[0].coefficient({G.domain().gen(0):G.degree()}))

        #monic centered model
        bool, L_id = field_in_database_NF(L)
        assert(bool)
        query['monic_centered.coeffs'] = [get_coefficients(g) for g in G]
        query['monic_centered.resultant'] = str(G.resultant())
        if L.degree() == 1:
            bad_primes = G.primes_of_bad_reduction()
        else:
            bad_primes = list(set([p.norm() for p in G.primes_of_bad_reduction()])) #remove duplicates
            bad_primes.sort()
        query['monic_centered.bad_primes'] = [int(p) for p in bad_primes]
        query['monic_centered.height'] = float(G.global_height())
        query['monic_centered.base_field_label'] = L_id
        #query['monic_centered.base_field.degree'] = int(L.degree())
        #models['original'].update({'base_field_emb': int(emb_index)})
        #conjugation to original model
        query['monic_centered.conjugation_from_original'] = [str(t) for r in M for t in r]
        #can these fields be different?
        query['monic_centered.conjugation_from_original_base_field_label'] = L_id

        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET monic_centered.coeffs = %(monic_centered.coeffs)s,
                monic_centered.resultant = %(monic_centered.resultant)s,
                monic_centered.bad_primes = %(monic_centered.bad_primes)s,
                monic_centered.base_field_label = %(monic_centered.base_field_label)s,
                monic_centered.conjugation_from_original = %(monic_centered.conjugation_from_original)s,
                monic_centered.conjugation_from_original_base_field_label = %(monic_centered.conjugation_from_original_base_field_label)s
            WHERE
                label = %(label)s
            """, query)


        log_file.write('monic centered model computed: ' + label + '\n')
        cancel_alarm()
        return True

    except AlarmInterrupt:
        log_file.write('monic centered model timeout: ' + str(timeout) + ':' + label + '\n')
    except:
        log_file.write('monic centered model failure: ' + label + '\n')
        raise
    cancel_alarm()
    return False

def add_chebyshev_model_NF(label, model_name='original', log_file=sys.stdout, timeout=30):
    """
    Determine if chebyshev and compute the chebyshev model.

    if 'is_polynomial' is not set, add_is_polynomial is called.

    if 'is_pcf' is not set, add_is_pcf is called.

    is_Chebshyev is set

    coeffs      varchar[],
    resultant   varchar,
    bad_primes  integer[],
    height      double precision,
    base_field_label varchar,
    conjugation_from_original varchar[],
    conjugation_from_original_base_field_label varchar(%s)

    """
    query={}
    query['label']=label
    log_file.write('starting chebyshev model for:' + label + '\n')
    my_cursor.execute("""SELECT is_polynomial FROM functions_dim_1_NF where label = %(label)s""",query)
    is_poly= my_cursor.fetchone()['is_polynomial']
    if is_poly is None:
        add_is_polynomial_NF(label, model_name=model_name, log_file=log_file, timeout=timeout)
        my_cursor.execute("""SELECT is_polynomial FROM functions_dim_1_NF where label = %(label)s""",query)
        is_poly= my_cursor.fetchone()['is_polynomial']
    if not is_poly:
        cancel_alarm()
        #not chebyshev if not a polynomial
        query['is_chebyshev']=False
        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET is_chebyshev = %(is_chebyshev)s
            WHERE
                label = %(label)s
            """, query)
        #TODO check rowcount
        log_file.write('chebyshev model done for:' + label + '\n')
        return True
    #check if chebyshev - Milnor
    my_cursor.execute("""SELECT is_pcf FROM functions_dim_1_NF where label = %(label)s""",query)
    is_pcf= my_cursor.fetchone()['is_pcf']
    if is_pcf is None:
        add_is_pcf(label, model_name=model_name, log_file=log_file, timeout=timeout)
        my_cursor.execute("""SELECT is_pcf FROM functions_dim_1_NF where label = %(label)s""",query)
        is_pcf= my_cursor.fetchone()['is_pcf']
    if not is_pcf:
        cancel_alarm()
        #not chebyshev if not pcf
        query['is_chebyshev']=False
        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET is_chebyshev = %(is_chebyshev)s
            WHERE
                label = %(label)s
            """, query)
        #TODO check rowcount
        log_file.write('chebyshev model done for:' + label + '\n')
        return True
    try:
        if timeout != 0:
            alarm(timeout)
        F = get_sage_func_NF(label, model_name=model_name, log_file=log_file)
        try:
            Fbar = F.change_ring(QQbar)
        except ValueError:
            Fbar = F.change_ring(F.base_ring().embeddings(QQbar)[0])
        d = F.degree()
        Pbar = Fbar.domain()
        crit, post_crit = get_post_critical(Fbar)
        count_crit = 0
        for Q in crit:
            if Q != Pbar([1,0]):
                count_crit += 1
        if count_crit == d-1:
            count = 0
            for Q in post_crit:
                if Q != Pbar([1,0]) and Q not in crit:
                    good_Q = Q
                    count += 1
            if count == 2:
                # and we need either [1,1] or both fixed
                m,n = good_Q.is_preperiodic(Fbar, return_period=True)
                if n == 1:
                    is_cheby = True
                else:
                    is_cheby = False
            else:
                is_cheby = False
        else:
            is_cheby = False
        query['is_chebyshev'] = is_cheby
        if not is_cheby:
            log_file.write('not chebyshev:' + label + '\n')
            my_cursor.execute("""UPDATE functions_dim_1_NF
            SET is_chebyshev = %(is_chebyshev)s
            WHERE
                label = %(label)s
            """, query)
            #TODO check rowcount
            return True
        #else compute a chebyshev model
        log_file.write('computing chebyshev model for:' + label + '\n')
        cheby_model = {}
        ch = F.domain().chebyshev_polynomial(d)
        conj_set = Fbar.conjugating_set(ch.change_ring(QQbar))
        K = ch.base_ring()
        bool, K_id = field_in_database_NF(K)
        assert(bool)

        query['is_chebyshev'] = True
        query['chebyshev_model.coeffs'] = [get_coefficients(g) for g in ch]
        query['chebyshev_model.resultant'] = str(ch.resultant())
        if K.degree() == 1:
            bad_primes = ch.primes_of_bad_reduction()
        else:
            bad_primes = list(set([p.norm() for p in ch.primes_of_bad_reduction()])) #remove duplicates
            bad_primes.sort()
        query['chebyshev_model.bad_primes'] = [int(p) for p in bad_primes]
        query['chebyshev_model.height'] = float(ch.global_height())
        query['chebyshev_model.base_field_label'] = K_id
        #query['chebyshev_model.base_field.degree'] = int(K.degree())
        #models['original'].update({'base_field_emb': int(emb_index)})
        #conjugation to original model
        N = ch.domain().dimension()
        M = conj_set[0]
        K, el, psi = number_field_elements_from_algebraics([t for r in M for t in r])
        L, phi = normalize_field_NF(K, log_file=log_file)
        M = matrix(L, N+1, N+1, [phi(t) for t in el])
        bool, L_id = field_in_database_NF(L)
        assert(bool)
        query['chebyshev_model.conjugation_from_original'] = [str(t) for r in M for t in r]
        query['chebyshev_model.conjugation_from_original_base_field_label'] = L_id

        #return query

        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET chebyshev_model.coeffs = %(chebyshev_model.coeffs)s,
                chebyshev_model.resultant = %(chebyshev_model.resultant)s,
                chebyshev_model.bad_primes = %(chebyshev_model.bad_primes)s,
                chebyshev_model.height = %(chebyshev_model.height)s,
                chebyshev_model.base_field_label = %(chebyshev_model.base_field_label)s,
                chebyshev_model.conjugation_from_original = %(chebyshev_model.conjugation_from_original)s,
                chebyshev_model.conjugation_from_original_base_field_label = %(chebyshev_model.conjugation_from_original_base_field_label)s,
                is_chebyshev = %(is_chebyshev)s
            WHERE
                label = %(label)s
            """, query)

        cancel_alarm()
        log_file.write('chebyshev model done for:' + label + '\n')
        #TODO check rowcount
        return True


    except AlarmInterrupt:
        log_file.write('chebyshev model timeout: ' + str(timeout) + ':' + label + '\n')
    except:
        log_file.write('chebyshev model failure: ' + label + '\n')
        raise
    cancel_alarm()
    return False

def add_newton_model_NF(label, model_name='original', log_file=sys.stdout, timeout=30):
    """
    Determine if newton and compute the newton model.
    See 1510.02271 for a possible newton citation

    'is_newton'
    coeffs      varchar[],
    resultant   varchar,
    bad_primes  integer[],
    height      double precision,
    base_field_label  varchar,
    polynomial_coeffs  varchar[],
    conjugation_from_original varchar[],
    conjugation_from_original_base_field_label varchar(%s)
    """
    query={}
    query['label']=label
    log_file.write('starting newton model for:' + label + '\n')
    #check if newton map
    try:
        F = get_sage_func_NF(label,model_name=model_name)
        N = F.domain().dimension()
        if timeout != 0:
            alarm(timeout)
        sigma_1 = F.sigma_invariants(1)
        d = ZZ(F.degree())
        newton_sigma = [d/(d-1)] + [0 for _ in range(d)] #almost newton
        if sigma_1 != newton_sigma:
            query['is_newton'] = False
            my_cursor.execute("""UPDATE functions_dim_1_NF
                SET is_newton = %(is_newton)s
                WHERE
                    label = %(label)s
                """, query)
            #TODO check rowcount
            cancel_alarm()
            log_file.write('newton model done for:' + label + '\n')
            return True
        #else is newton
        try:
            Fbar = F.change_ring(QQbar)
        except ValueError:
            Fbar = F.change_ring(F.base_ring().embeddings(QQbar)[0])
        Pbar = Fbar.domain()
        fixed = Fbar.periodic_points(1)
        for Q in fixed:
            if Fbar.multiplier(Q,1) != 0:
                inf = Q
                break
        if inf != Pbar([1,0]):
            #need to move inf to infinity
            fixed.remove(inf)
            source = [inf] + fixed[:2]
            target = [Pbar([1,0]), Pbar([0,1]), Pbar([1,1])]
            M = Pbar.point_transformation_matrix(source, target)
            M = M.inverse()
            newton = Fbar.conjugate(M)
            K, el, psi = number_field_elements_from_algebraics([t for r in M for t in r])
            L, phi = normalize_field_NF(K, log_file=log_file)
            M = matrix(L, N+1, N+1, [phi(t) for t in el])

            newton = newton._number_field_from_algebraics()
            #normalize base field
            newton, phi = normalize_function_NF(newton, log_file)
            #fix variable names
            PN = ProjectiveSpace(newton.base_ring(), Pbar.dimension(), Pbar.variable_names())
            RN = PN.coordinate_ring()
            newton = DynamicalSystem([RN(newt) for newt in newton], domain=PN)
        else:
            newton = F
            L = F.base_ring()
            M = matrix(QQ,2,2,[1,0,0,1])
        bool, L_id = field_in_database_NF(L)
        assert(bool)

        N_aff = newton.dehomogenize(1)
        z = N_aff.domain().gen(0)
        Npoly = (z-N_aff[0]).numerator()
        assert(Npoly.derivative(z) == (z-N_aff[0]).denominator()), "not actually newton"
        query['is_newton'] = True

        query['newton_model.coeffs'] = [get_coefficients(g) for g in newton]
        query['newton_model.resultant'] = str(newton.resultant())
        if L.degree() == 1:
            bad_primes = newton.primes_of_bad_reduction()
        else:
            bad_primes = list(set([p.norm() for p in newton.primes_of_bad_reduction()])) #remove duplicates
            bad_primes.sort()
        query['newton_model.bad_primes'] = [int(p) for p in bad_primes]
        query['newton_model.height'] = float(newton.global_height())
        query['newton_model.base_field_label'] = L_id
        #query['newton_model.base_field.degree'] = int(L.degree())
        #models['original'].update({'base_field_emb': int(emb_index)})
        #conjugation to original model
        query['newton_model.conjugation_from_original'] = [str(t) for r in M for t in r]
        query['newton_model.conjugation_from_original_base_field_label'] = L_id
        C = []
        z = Npoly.parent().gen(0)
        for i in range(0,Npoly.degree()+1):
            C.append(str(Npoly.coefficient({z:i})))
        query['newton_model.polynomial_coeffs'] = C

        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET newton_model.coeffs = %(newton_model.coeffs)s,
                newton_model.resultant = %(newton_model.resultant)s,
                newton_model.bad_primes = %(newton_model.bad_primes)s,
                newton_model.height = %(newton_model.height)s,
                newton_model.base_field_label = %(newton_model.base_field_label)s,
                newton_model.conjugation_from_original = %(newton_model.conjugation_from_original)s,
                newton_model.conjugation_from_original_base_field_label = %(newton_model.conjugation_from_original_base_field_label)s,
                newton_model.polynomial_coeffs = %(newton_model.polynomial_coeffs)s,
                is_newton = %(is_newton)s
            WHERE
                label = %(label)s
            """, query)

        cancel_alarm()
        log_file.write('newton model done for:' + label + '\n')
        #TODO check rowcount
        return True

    except AlarmInterrupt:
        log_file.write('newton model timeout: ' + str(timeout) + ':' + label + '\n')
    except:
        log_file.write('newton model failure: ' + label + '\n')
        raise
    cancel_alarm()
    return False


def add_is_lattes_NF(label, model_name='original', log_file=sys.stdout, timeout=30):
    """
    Determine if lattes. Updates is_pcf if not known

    'is_lattes'
    """
    log_file.write('starting is lattes for:' + label + '\n')
    query={}
    query['label']=label
    # must be pcf
    my_cursor.execute("""SELECT is_pcf FROM functions_dim_1_NF where label = %(label)s""",query)
    is_pcf= my_cursor.fetchone()['is_pcf']
    if is_pcf is None:
        add_is_pcf(label, model_name=model_name, log_file=log_file, timeout=timeout)
        my_cursor.execute("""SELECT is_pcf FROM functions_dim_1_NF where label = %(label)s""",query)
        is_pcf= my_cursor.fetchone()['is_pcf']
    if not is_pcf:
        cancel_alarm()
        #not lattes if not pcf
        query['is_lattes']=False
        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET is_lattes = %(is_lattes)s
            WHERE
                label = %(label)s
            """, query)
        #TODO check rowcount
        log_file.write('lattes done for:' + label + '\n')
        return True

    #check if lattes map
    try:
        if timeout != 0:
            alarm(timeout)

        F = get_sage_func_NF(label,model_name=model_name, log_file=log_file)
        d = ZZ(F.degree())
        try:
            Fbar = F.change_ring(QQbar)
        except ValueError:
            Fbar = F.change_ring(F.base_ring().embeddings(QQbar)[0])
        Pbar = Fbar.domain()

        crit, post_crit = get_post_critical(Fbar)
        if (len(crit) == 2*d - 2) and \
            (len(set(post_crit).difference(set(crit))) == 4):
            is_lattes = True
            #TODO get curve
            """
            try:
                if timeout != 0:
                    alarm(timeout)
                P2.<z,w> = ProjectiveSpace(QQbar,1)
                T = 1
                for t in post_crit:
                    T = T*(z-t[0])
                C = Curve(w^2-T)
                Q = C(post_crit[0][0],0)
                E = EllipticCurve_from_plane_curve(C,Q)
                mult_1 = F.multiplier_spectra(1)
                for t in mult_1:#multipliers are m, -m, m^2
                    if t^2 in mult_1:
                        m = t
                        break
                m = max([-m,m])
                f['lattes_info'] = {'curve': str(E.defining_polynomial()), 'LMFDB_label' : 'xxxx', 'm' : str(m)}
            except AlarmInterrupt:
                log_file.write('timeout lattes: ' + str(timeout) + ':' + models['original']['polys']['val'] + '\n')
            """
        else:
           is_lattes = False

        query['is_lattes'] = is_lattes
        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET is_lattes = %(is_lattes)s
            WHERE
                label = %(label)s
            """, query)
        log_file.write('is lattes complete for ' + label + '\n')
        cancel_alarm()
        #TODO check rowcount
        return True

    except AlarmInterrupt:
        log_file.write('is lattes timeout: ' + str(timeout) + ':' + label + '\n')
    except:
        log_file.write('is lattes failure: ' + label + '\n')
        raise
    cancel_alarm()
    return False

def add_citations_NF(label, citations, log_file=sys.stdout):
    """
        Add the id of the citations for this functions
    """
    #make sure they are ints
    citations = [int(t) for t in citations]

    my_cursor.execute("""SELECT
        citations
         FROM functions_dim_1_NF
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
    my_cursor.execute("""UPDATE functions_dim_1_NF
        SET citations = %s
        WHERE
            label = %s
        """, [new_cites, label])

def add_function_all_NF(F, log_file=sys.stdout):
    """
    add all entries for one dynamical system
    """
    label=add_function_NF(F,log_file=log_file)
    add_citations_NF(label, citations, log_file=log_file)
    add_sigma_inv_NF(label,'original',2,3,log_file=log_file)
    add_is_pcf(label,'original',log_file=log_file)
    add_critical_portrait(label,'original',log_file=log_file)
    add_automorphism_group_NF(label,'original',log_file=log_file)
    add_rational_preperiodic_points_NF(label,log_file=log_file)
    if label not in ['1.2.c8ddd2a7.1', '1.2.611af5d0.1', '1.2.7cb63904.1',\
        '1.2.0d94343e.1','1.2.beb83dd2.1','1.2.dd6fd9ae.1','1.2.4994c36e.1',\
        '1.2.98d76bdd.1']:
        add_reduced_model_NF(label,log_file=log_file)
    add_is_polynomial_NF(label,log_file=log_file)
    add_monic_centered_model_NF(label,log_file=log_file)
    add_chebyshev_model_NF(label,log_file=log_file)
    add_newton_model_NF(label,log_file=log_file)
    add_is_lattes_NF(label,log_file=log_file)
    choose_display_model(label,log_file=log_file)

    return label

