"""
Functionality for working with number fields

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


import sys
from sage.categories.number_fields import NumberFields
from sage.libs.pari import pari
from sage.rings.cc import CC
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from six.moves.urllib.request import urlopen
import json


def normalize_field_NF(K, emb=None, log_file=sys.stdout):
    """
    put the base field in normalized form.
    Return the normalized field and an embedding
    into it from the given field. If we are given a number field
    with no embedding, then one is created at random.
    """
    if K in NumberFields():
        if K != QQ:
            def_poly = K.defining_polynomial()
            def_poly_new, gen = pari(def_poly).polredabs(1)
            def_poly_new = def_poly_new.sage({def_poly.parent().variable_names()[0]:def_poly.parent().gen(0)})
            if def_poly_new != def_poly:
                gen = gen.lift().sage({def_poly.parent().variable_names()[0]:def_poly.parent().gen(0)})
                if emb is None:
                    if K.coerce_embedding() is None:
                        L = NumberField(def_poly_new, 'a', embedding=CC(def_poly_new.roots(ring=CC)[0][0]))
                    else:
                        L = NumberField(def_poly_new, 'a', embedding=gen(K.coerce_embedding()(K.gen())))
                else:
                    L = NumberField(def_poly_new, 'a', embedding=gen(emb(K.gen())))
                to_L = K.hom([gen], L)
                log_file.write('Normalized field:' + str(K) + ' to ' + str(L) + '\n')
                return L, to_L
            if K.variable_name() != 'a':
                L = K.change_names('a')
                log_file.write('Normalized field variable:' + str(K) + ' to ' + str(L) + '\n')
                return L, L.structure()[1]
        return K, K.automorphisms()[0].parent().identity()
    #other field
        raise NotImplementedError('Can only normalize number fields with this function')


def check_field_normalized_NF(K, log_file=sys.stdout):
    """
    check that base field is in normalize form.

    For number fields this is a monic polynomial.
    The generator is named 'a'.
    """
    if K is QQ:
        return True

    if K in NumberFields():
        if not K.is_absolute():
            return False
        def_poly = K.defining_polynomial()
        if K.variable_name() != 'a':
            return False
        def_poly_new = pari(def_poly).polredabs()
        if def_poly == pari(def_poly).polredabs():
            return True
        return False

    raise ValueError('field not a number field')


def lmfdb_field_label_NF(K, log_file = sys.stdout): #label finding
    #assumes the field is nomalized
    if not check_field_normalized_NF(K):
        raise ValueError('field not normalized')
    poly = K.defining_polynomial() 
    z = poly.parent().gen(0)
    C = poly.coefficients(sparse=False)
    url = 'https://beta.lmfdb.org/api/nf_fields/?coeffs={'
    for i in range(len(C)-1):
        url += str(C[i]) + ','
    url += str(C[-1]) + '}&_format=json&_fields=label&_delim=;'
    page = urlopen(url)
    dat = str(page.read().decode('utf-8'))
    dat = json.loads(dat)['data']
    if dat != []:
        label = dat[0]['label']
        return True, label
    else:
        # can't find the field
        log_file.write('Field not found in LMFDB: ' + str(K) + '\n')
        label = 0
        return False, label


def get_sage_field_NF(label): # get field from db, return as sage object

    url = 'https://beta.lmfdb.org/api/nf_fields/?label=' + label + '&_format=json&_fields=coeffs&_delim=;'
    page = urlopen(url)
    dat = str(page.read().decode('utf-8'))
    dat = json.loads(dat)['data']
    if dat == []:
        log_file.write('Field not found in LMFDB: ' + str(label) + '\n')
        return 0
    C = dat[0]['coeffs']
    if C == [0,1]:
        return QQ
    R=PolynomialRing(QQ, 'z')
    z = R.gen(0)
    poly = 0
    for i in range(len(C)):
        poly += z**i*C[i]
    K = NumberField(poly, 'a')
    return K
