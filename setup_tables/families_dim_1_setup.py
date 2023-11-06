"""
Setup PostgreSQL table for families of functions in dimension 1.

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
    DROP TABLE IF EXISTS families_dim_1_NF
""")


############################
# create custom types

# create new types


######################################
# Create Functions table Schema

my_cursor.execute("""
CREATE TABLE families_dim_1_NF (
    label varchar(%s) PRIMARY KEY,
    degree integer,
    num_parameters integer,
    model_coeffs varchar[],
    model_resultant varchar,
    base_field_label varchar(%s),
    base_field_degree integer,
    sigma_invariants sigma_invariants_type,
    citations integer[],
    is_polynomial boolean,
    num_critical_points integer,
    automorphism_group_cardinality integer
    )
""",[function_label_length, field_label_length])


my_session.commit()
