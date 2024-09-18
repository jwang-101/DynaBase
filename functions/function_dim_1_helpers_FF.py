"""
Functions to help working with functions in dimension 1
defined over finite fields.

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
import sys

#length of SHAKE-256S hash for function label
digest_length = int(4)


###########################################
# General Functions
###########################################


from functions.function_dim_1_helpers_generic import get_coefficients
from functions.function_dim_1_helpers_generic import get_post_critical
from functions.function_dim_1_helpers_generic import choose_display_model
from functions.function_dim_1_helpers_generic import graph_to_array
from functions.function_dim_1_helpers_generic import array_to_graph


##############################################################
#  Functionality for working with functions over finite fields
##############################################################

def normalize_function_FF(F, log_file=sys.stdout):
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
        K, phi = normalize_field_FF(base_field, log_file=log_file)
        if not K is base_field:
            #extend to polyring
            phi = F.base_ring().hom(phi, F.base_ring().change_ring(phi.codomain()))
            F = F.change_ring(phi)
            log_file.write('Normalized base field:' + str(F.base_ring()) + ' to ' + str(K) + '\n')
    else:
        K, phi = normalize_field_FF(base_field, log_file=log_file)
        if not K is base_field:
            F = F.change_ring(phi)
            log_file.write('Normalized base field:' + str(base_field) + ' to ' + str(K) + '\n')

    return F, phi

    #########these need to be updated#################

def get_sage_func_FF(model_name, label=None, f=None, log_file=sys.stdout):
    """
    Given a database entry f or label and a model name, return the
    sage dynamical system.

    todo: base field is family, i.e., polynomial ring or function field
    """
    if not label is None:
        f = functions.find_one({'label':label})
    if f is None:
        raise ValueError('no function specified')

    base_field = get_sage_field(f['models'][model_name]['base_field_label'], f['base_field_type'],\
                emb_index=f['models'][model_name]['base_field_emb'])
    num_parameters = f['models'][model_name]['num_parameters']
    if num_parameters != 0:
        R = PolynomialRing(base_field, num_parameters, 't')
        P = ProjectiveSpace(f['N'], R, f['models'][model_name]['polys']['vars'])
        var_dict = P.coordinate_ring().gens_dict()
        var_dict.update(R.gens_dict())
    else:
        P = ProjectiveSpace(f['N'], base_field, f['models'][model_name]['polys']['vars'])
        var_dict = P.coordinate_ring().gens_dict()
    var_dict.update(base_field.gens_dict())
    polys = sage_eval(f['models'][model_name]['polys']['val'], var_dict)
    return DynamicalSystem(polys, domain=P)

def check_conjugates_FF(F,G, normalize_base=False, field_type=None, num_parameters=None):
    """
    F,G are two sage models with same degree, dimension

    returns 0 for different conj class
    1 for rational twists
    2 for rationally conjuate

    """
    if num_parameters > 0:
        raise NotImplementedError('cannot check conjugates for families')
        log_file.write('cannot check conjugates for families: ' + str(list(F)) + ' : ' + str(list(G)) + '\n')
    if field_type == 'p-adic field':
        raise NotImplementedError('cannot check conjugates for p-adic fields')
        log_file.write('cannot check conjugates for p-adic fields: ' + str(list(F)) + ' : ' + str(list(G)) + '\n')
    if normalize_base:
        F, phiF = normalize_function(F)
        G, phiG = normalize_function(G)
    Kf = F.base_ring()
    Kg = G.base_ring()

    if Kf.degree() > 1:
        raise NotImplementedError('cannot check conjugates for non-prime finite fields')
        log_file.write('cannot check conjugates for non-prime finite fields: ' + str(list(F)) + ' : ' + str(list(G)) + '\n')
    KFbar = Kf.algebraic_closure()
    Fbar = F.change_ring(KFbar)
    KGbar = Kg.algebraic_closure()
    Gbar = G.change_ring(KFbar)
    CS = Fbar.conjugating_set(Gbar)

    if len(CS) == 0:
        return 0
    Kf = F.base_ring()
    Kg = G.base_ring()
    if Kf != Kg:
        return 1 #either conjugate or not, but not rational twists since defined over different fields
    #so we assume they have the same base field so we are looking for twists
    for m in CS:
        field_list = [t.as_finite_field_element()[0] for r in m for t in r]
        normalized_field_list = [normalize_field(L, log_file=log_file)[0] for L in field_list]
        if all([K.degree() <= Kf.degree() for K in normalized_field_list]):
            return 1
    return 2
                                                        # unsure if this max_sigma counts as a sigma3 thingy
def conj_in_database_FF(F, sigma_1=None, conj_fns=None, max_sigma=3, compare_model='original', log_file=sys.stdout, timeout=30):
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
    if isinstance(base_field, (PolynomialRing_general, MPolynomialRing_base)) or base_ring in FunctionFields():
        num_parameters = int(base_field.ngens())
        base_field = F.base_ring().base_ring()
    else:
        num_parameters = int(0)
        base_field = F.base_ring()
    field_label, emb_index, field_id = get_field_label(base_field, log_file=log_file)

    if conj_fns is None:
        if sigma_1 is None:
            sigma_1 = str(F.sigma_invariants(1))
        query = {'degree':int(F.degree()), 'N':int(F.domain().dimension_relative()), 'sigma_inv.1': sigma_1}
        query['base_field_label'] = field_label
        query['num_parameters'] = num_parameters
        conj_fns = functions.find(query)
    else:
        query = {}
    polys = {'vars':list(F.domain().variable_names()), 'val':str(list(F)), 'latex':latex(list(F))}

    bool, g = model_in_database(F, conj_fns=conj_fns, log_file=log_file)
    if bool:
        return 1,g
    #now check for twists/conjugates
    conj_fns.rewind()
    if query.has_key('sigma_inv'):
        s = len(query['sigma_inv'].keys)
    else:
        s = 1
        query = {'degree':int(F.degree()), 'N':int(F.domain().dimension_relative())}
        query['base_field_label'] = field_label
        query['num_parameters'] = num_parameters
        query['sigma_inv'] = {}
    while s <= max_sigma and conj_fns.count() != 0:
        #need to be careful about querying too high
        if conj_fns.count() != 0:
            max_s = max([len(g['sigma_inv'].keys()) for g in conj_fns])
            conj_fns.rewind()
        if s <= max_s:
            query['sigma_inv'].update({str(s): str(F.sigma_invariants(s))})
        else:
            log_file.write('trying too many sigma invariants : ' + str(s) + ' > ' + str(max_s) + '\n')
        conj_fns = functions.find(query)
        s += 1
    if isinstance(base_field, (PolynomialRing_general, MPolynomialRing_base)) or base_ring in FunctionFields():
        base_field = base_field.base_ring()
    K = fields.find_one({'label':field_label})
    try:
        if timeout != 0:
            alarm(timeout)
        for g in conj_fns:
            twist_val = check_conjugates(F, get_sage_func(compare_model, f=g),\
                                field_type=K['type'], num_parameters=num_parameters)
            if twist_val == 1:
                log_file.write('already there: conjugate found in db: ' + polys['val'] + ' as ' + g['models']['original']['polys']['val'] + '\n')
                cancel_alarm()
                return 1, g
            elif twist_val == 2:
                if g.has_key('twists'):
                    twist_list = g['twists'] + [g['label']]
                    twist_list.sort()
                else:
                    twist_list = [g['label']]
                log_file.write('has twist in db: ' + polys['val'] + ' as ' + str(twist_list) + '\n')
                cancel_alarm()
                return 2, twist_list
        #not in database
        cancel_alarm()
        return 0, []

    except AlarmInterrupt:
        log_file.write('timeout: func_in_db: ' + str(timeout) + ':' + polys['val'] + '\n')
        raise

def add_function_FF(F, citations=[], family=[], keywords=[], action='add', bool_add_field=False, log_file=sys.stdout, timeout=30):
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

    #todo need to update conjugate check:
    #This checks if here are any twists or congjuates in the database
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

    log_file.write('Searching for functions: ' + label + ': ')
    query = {'degree':f['degree'], 'sigma_invariants.one': f['sigma_invariants.one']}
    my_cursor.execute("""SELECT * FROM functions_dim_1_NF
        WHERE degree=%(degree)s AND
            (sigma_invariants).one=%(sigma_invariants.one)s::varchar[]""",query)
    conj_fns = my_cursor.fetchall()
    m = len(conj_fns)
    found = 0 #hard coded until check updated
    #found, L = conj_in_database_NF(F, conj_fns=conj_fns, log_file=log_file, timeout=timeout)
    L={} #hard coded until check updated

    # TODO: need to update found == 1 section
    if found == 1:
        #function or rational conjugate already there
        if action == 'add':
            log_file.write('function already known for : ' + L['label'] + '\n')
            return L[label]
        elif action == 'replace':
            #delete old function
            delete_one_function_NF(L['label'], log_file=log_file)
            log_file.write('replacing function : ' + L['label'] + '\n')
            #adjust the count for deletion
            #what do we do about deleting a label in the middle?
            m = m-1
        else:
            raise ValueError('action must be add/replace')
   #otherwise we'll add the function
   #assume none are conjugate if we get to here
    f['label'] = label + str(m+1)
    #TODO: update citation code
    cites = []
    for cite in citations:
        cites.append(bibliography[cite])
    if cites == []:
        f['citations'] = None
    else:
        f['citations'] = cites
    if family == []:
        f['family'] = None
    else:
        f['family'] = family

    #original model
    f['original_model.coeffs'] = [get_coefficients(g) for g in F]
    f['original_model.resultant'] = str(F.resultant())
    if F.base_ring().degree() == 1:
        bad_primes = F.primes_of_bad_reduction()
    else:
        bad_primes = list(set([p.norm() for p in F.primes_of_bad_reduction()])) #remove duplicates
        bad_primes.sort()
    f['original_model.bad_primes'] = bad_primes
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
        (label, degree, base_field_label,
         base_field_degree, sigma_invariants.one, citations, family,
         original_model.coeffs, original_model.resultant, original_model.bad_primes,
         original_model.base_field.label, original_model.base_field.degree,
         original_model.conjugation_from_original, original_model.conjugation_from_original_base_field_label)
        VALUES
        (%(label)s, %(degree)s, %(base_field_label)s,
         %(base_field_degree)s, %(sigma_invariants.one)s, %(citations)s, %(family)s,
         %(original_model.coeffs)s, %(original_model.resultant)s, %(original_model.bad_primes)s,
         %(original_model.base_field_label)s,
         %(original_model.conjugation_from_original)s, %(original_model.conjugation_from_original_base_field_label)s)
        RETURNING label """,f)
    F_id = my_cursor.fetchone()[0]
    log_file.write('inserted: ' + str(list(F)) + ' as ' + str(F_id) + '\n')
    # TODO: update twists check
    if found == 2:
        #There is one or more twists in database
        f['twists'] = L
        log_file.write('twist list: ' + str(L) + ' for ' + f['label'] + '\n')
        functions.update_one({'_id':ins_id},{'$set' : {'twists' : L}})
        for twist_label in L:
            h = functions.find_one({'label':twist_label})
            twist_list = copy(L)
            twist_list.remove(twist_label)
            twist_list.append(f['label'])
            twist_list.sort()
            functions.update_one({'_id':h['_id']},{'$set' : {'twists' : twist_list}})
            log_file.write('updating twist list for: ' + h['label'] + '\n')
    return F_id



def add_sigma_inv_FF(label, model_name='original', start=2, end=3, log_file=sys.stdout, timeout=30):
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
            # # i think im supposed to delete the below code?
            # elif k == 3:
            #     sigma.update({'three':[str(t) for t in F.sigma_invariants(k)]})
            #     my_cursor.execute("""UPDATE functions_dim_1_NF
            #         SET sigma_invariants.three = %(three)s
            #         WHERE label=%(label)s
            #         """,sigma)
            #     if my_cursor.rowcount == 0:
            #         log_file.write('function ' + label + 'not found \n')
            #         raise ValueError("function not found to update")
            #     else:
            #         log_file.write('updated ' + str(my_cursor.rowcount) + ' functions for sigma.' + str(k) + '\n')
            else:
                raise NotImplementedError("only upto k=3") #should i edit this only up to k = 2?

        cancel_alarm()
        return True
    except (AlarmInterrupt, KeyboardInterrupt):
        log_file.write('timeout: sigma inv for:' + str(timeout) + ':' + label + '\n')
    except:
        log_file.write('failure: sigma inv for:' + label + '\n')
        raise

    cancel_alarm()
    return False


def add_automorphism_group_FF(label, model_name='original', log_file=sys.stdout, timeout=30):
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

def add_rational_preperiodic_points_FF(label, model_name='original', field_label=None, log_file=sys.stdout, timeout=30):
    """
    Find the rational preperiodic points and add

    number_rational_preperiodic integer,
    rational_periodic_cycles integer[],
    rational_preperiodic_components integer[],
    rational_preriodic_graph_id varchar,
    rational_periodic_points varchar[][2],
    """
    if timeout != 0:
        alarm(timeout)

    log_file.write('starting rational preperiodic points for:' + label + '\n')
    try:
        query={}
        query['label']=label
        F = get_sage_func_NF(label, model_name=model_name, log_file=log_file)
        if field_label is None:
            K = F.base_ring()
        else:
            #TODO move preperiodic data to it's own table to account for other fields
            #for the same function
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

        query['number_rational_preperiodic'] = len(preper.vertices())
        query['rational_periodic_cycles'] = [int(len(t)-1) for t in preper.all_simple_cycles()]
        query['rational_preperiodic_components'] = preper.connected_components_sizes()
        #one per component
        query['rational_periodic_points'] = [[str(t) for t in T[0]] for T in g.all_simple_cycles()]
        my_cursor.execute("""UPDATE functions_dim_1_NF
            SET number_rational_preperiodic = %(number_rational_preperiodic)s,
                rational_periodic_cycles = %(rational_periodic_cycles)s,
                rational_preperiodic_components = %(rational_preperiodic_components)s,
                rational_periodic_points = %(rational_periodic_points)s
            WHERE label=%(label)s
            """,query)
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

def add_is_polynomial_FF(label, model_name='original', log_file=sys.stdout, timeout=30):
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

def add_monic_centered_model_FF(label, model_name='original', log_file=sys.stdout, timeout=30):
    """
    Compute the monic centered model.

    If it hasn't already been checked to be, add_is_polynomial_NF is called.

    coeffs      varchar[],
    resultant   varchar,
    bad_primes  integer[],
    height      double precision,
    base_field  varchar,
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
        add_is_polynomial(label, model_name=model_name, log_file=log_file, timeout=timeout)
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
            #elif f['base_field_type'] == 'finite field':
            #g,M,phi = F.normal_form(return_conjugation=True)
            ##base field may have changed so normalize again
            #g, phi = normalize_function(g, log_file)
            #L = g.base_ring()
            #M = matrix(L, f['N']+1, f['N']+1, [phi(t) for r in M for t in r])

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
        query['monic_centered.bad_primes'] = bad_primes
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


def add_chebyshev_model_FF(label, model_name='original', log_file=sys.stdout, timeout=30):
    """
    Determine if chebyshev and compute the chebyshev model.

    if 'is_polynomial' is not set, add_is_polynomial is called.

    if 'is_pcf' is not set, add_is_pcf is called.

    is_Chebshyev is set

    coeffs      varchar[],
    resultant   varchar,
    bad_primes  integer[],
    height      double precision,
    base_field  varchar,
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
        is_poly= my_cursor.fetchone()['is_pcf']
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
                chebyshev_model.base_field_label = %(chebyshev_model.base_field_label)s,
                chebyshev_model.conjugation_from_original = %(chebyshev_model.conjugation_from_original)s,
                chebyshev_model.conjugation_from_original_base_field_label = %(chebyshev_model.conjugation_from_original_base_field_label)s
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
