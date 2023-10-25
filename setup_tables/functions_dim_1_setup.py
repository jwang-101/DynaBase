"""
Set up the PostgreSQL tables for functions in dimension 1
defined over number fields and finite fields.

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
    DROP TABLE IF EXISTS functions_dim_1_NF
""")

my_cursor.execute("""
    DROP TABLE IF EXISTS functions_dim_1_FF
""")

############################
# create custom types

my_cursor.execute("""
DROP TYPE IF EXISTS model_type CASCADE;
""")
my_cursor.execute("""
DROP TYPE IF EXISTS newton_model_type CASCADE;
""")
my_cursor.execute("""
DROP TYPE IF EXISTS sigma_invariants_type CASCADE;
""")
my_cursor.execute("""
DROP TYPE IF EXISTS base_field_type CASCADE;
""")


my_cursor.execute("""
DROP TYPE IF EXISTS display_model_type CASCADE;
""")


#create new types
my_cursor.execute("""
CREATE TYPE display_model_type AS ENUM ('original', 'reduced', 'monic centered', 'chebyshev', 'newton')
""")


my_cursor.execute("""
CREATE TYPE base_field_type AS (
    label   varchar(%s),
    degree  integer
  )
  """,
  [field_label_length])


my_cursor.execute("""
CREATE TYPE sigma_invariants_type AS (
    one      varchar[],
    two      varchar[],
    three    varchar[]
  )
""")

my_cursor.execute("""
CREATE TYPE model_type AS (
    coeffs      varchar[],
    resultant   varchar,
    bad_primes  integer[],
    height      double precision,
    base_field_label  varchar(%s),
    conjugation_from_original varchar[],
    conjugation_from_original_base_field_label varchar(%s)
  )""",
  [field_label_length,field_label_length])

my_cursor.execute("""
CREATE TYPE newton_model_type AS (
    coeffs      varchar[],
    resultant   varchar,
    bad_primes  integer[],
    height      double precision,
    base_field_label  varchar(%s),
    polynomial_coeffs  varchar[],
    conjugation_from_original varchar[],
    conjugation_from_original_base_field_label varchar(%s)
  )""",
  [field_label_length,field_label_length])


######################################
# Create Functions table Schema

#Need Newton polynomial for newton model
#rational preperiodic points should be one point per cycle
#  should the data type be a list of points or a list of strings?
#Should the rational preperiodic point information be put into a separate table
#  to make multiple entries per function more efficient?
#Can we find the elliptic curve of Lattes?

my_cursor.execute("""
CREATE TABLE functions_dim_1_NF (
    label varchar(%s) PRIMARY KEY,
    degree integer,
    base_field_label varchar(%s),
    base_field_degree integer,
    sigma_invariants sigma_invariants_type,
    citations integer[],
    family varchar[],
    original_model model_type,
    monic_centered model_type,
    chebyshev_model model_type,
    reduced_model model_type,
    newton_model newton_model_type,
    display_model display_model_type,
    is_polynomial boolean,
    is_chebyshev boolean,
    is_newton boolean,
    is_lattes boolean,
    is_pcf boolean,
    automorphism_group_cardinality integer,
    rational_twists varchar[],
    number_rational_preperiodic integer,
    rational_periodic_cycles integer[],
    rational_preperiodic_components integer[],
    rational_preriodic_graph_id varchar,
    rational_periodic_points varchar[][2],
    critical_portrait_cardinality integer,
    post_critical_cardinality integer,
    critical_portrait_components integer[],
    critical_portrait_structure integer[][2],
    critical_portrait_graph_id varchar
  )
""",[function_label_length, field_label_length])


my_cursor.execute("""
CREATE TABLE functions_dim_1_FF (
    label varchar(%s) PRIMARY KEY,
    degree integer,
    base_field_label varchar(%s),
    base_field_degree integer,
    sigma_invariants sigma_invariants_type,
    citations integer[],
    family varchar[],
    original_model model_type,
    monic_centered model_type,
    chebyshev_model model_type,
    reduced_model model_type,
    newton_model newton_model_type,
    display_model display_model_type,
    is_polynomial boolean,
    is_chebyshev boolean,
    is_newton boolean,
    is_lattes boolean,
    is_permutation boolean,
    automorphism_group_cardinality integer,
    rational_twists varchar[],
    rational_periodic_cycles integer[],
    rational_preperiodic_components integer[],
    avg_tail_length real
  )
""",[function_label_length,field_label_length])


my_session.commit()
