"""
Functionality for working with finite fields

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

def normalize_field_FF(K, log_file=sys.stdout):
    """
    put the base field in normalized form.
    Return the normalized field and an embedding
    into it from the given field. If we are given a number field
    with no embedding, then one is created at random.
    """
    if K in FiniteFields():
        CP = ConwayPolynomials()
        MR = K.modulus().parent()
        if MR(K.modulus()) != MR(CP.polynomial(K.characteristic(), K.degree())):
            raise NotImplementedError('must be conway modulus')
        if K.variable_name() != 'a':
            L = GF(K.cardinality(), 'a', modulus='conway')
            log_file.write('Normalized field variable:' + str(K) + ' to ' + str(L) + '\n')
            return L, K.hom([L.gen(0)], L)
        return K, K.Hom(K).identity()
    else: #other field
        raise NotImplementedError('Can only normalize finite fields with this function')
        

def check_field_normalized_FF(K, log_file=sys.stdout):
    """
    check that base field is in normalize form.

    For finite fields the conway polynomial is used
    """
    if K in FiniteFields():
        if K.is_prime_field():
            return True
        if K.variable_name() != 'a':
            return False
        CP = ConwayPolynomials()
        MR = K.modulus().parent()
        try:
            if K.modulus() == MR(CP.polynomial(K.characteristic(),K.degree())):
                return True
        except RuntimeError: #conway poly not in database
            return False
        return False

    raise ValueError('field not finite')


def field_in_database_FF(K, log_file=sys.stdout):
    if not check_field_normalized_FF(K, log_file=log_file):
        raise ValueError("field not normalized")
    D = {}
    D['modulus_coeffs'] = [int(t) for t in K.modulus().coefficients(sparse=False)]
    my_cursor.execute("""SELECT label FROM finite_fields WHERE modulus_coeffs=%(modulus_coeffs)s""",D)
    log_file.write('Field not Found: ' + str(K) + '\n')
    label = my_cursor.fetchone()
    if label is None:
        return False, 0
    return True, label['label']
    
def add_field_FF(K, normalize=True, log_file=sys.stdout):
    """
    Add the field K to the finite_fields table.

    finite field: p.deg

    TODO::this needs timeout
    """
    log_file.write('Adding field: ' + str(K) + ':')
    if not check_field_normalized_FF(K, log_file=log_file):
        if not normalize:
            raise ValueError('field not normalized')
        else:
            K, phi = normalize_field_FF(K, log_file=log_file)
    bool, K_id = field_in_database_FF(K, log_file=log_file)
    if bool:
        log_file.write('already in db: ' + str(K) + '\n')
        #already in database
        return -1
    
    if K in FiniteFields():
        label = str(K.characteristic()) + '.' + str(K.degree())
    else:
        raise NotImplementedError("only finite fields")

    log_file.write('label computed: ' + label + '\n')

    F = {}
    F['label'] = label
    F['modulus_coeffs'] = [int(t) for t in K.modulus().coefficients(sparse=False)]
    F['characteristic'] = int(K.characteristic())
    F['cardinality'] = int(K.cardinality())

    # add to finite fields
    my_cursor.execute("""INSERT INTO finite_fields VALUES
        (%(label)s, %(modulus_coeffs)s, %(characteristic)s,%(cardinality)s)
        RETURNING label """,F)
    K_id = my_cursor.fetchone()['label']
    log_file.write('field inserted as: ' + str(K_id) + '\n')
    return K_id


def delete_field_FF(K):
    bool, id = field_in_database_FF(K)
    if not bool:
        raise ValueError("field not in database")
    val = my_cursor.execute("""DELETE FROM finite_fields WHERE label = K_label""")
    if val==0:
        raise ValueError("Failed to find field" + str(K))
    elif val>1:
        raise ValueError("Deleted more than 1 row") #should never occur
    return val

def get_sage_field_FF(field_label):
    """
    given a field label, get the field from the database and return the
    field as sage object
    """
    D={'field_label':field_label}
    my_cursor.execute("""SELECT * FROM finite_fields WHERE label=%(field_label)s""",D)
    F = my_cursor.fetchone()
    R = PolynomialRing(QQ,'X')
    def_poly = R(F['defn_poly_coeffs'])
    base_field = GF(F['characteristic']**F['degree'], 'a', modulus=def_poly)
    return base_field
