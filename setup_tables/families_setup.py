"""
Setup PostgreSQL table for families of functions.

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
    DROP TABLE IF EXISTS families_dim_1
""")


############################
# create custom types

# create new types


######################################
# Create Functions table Schema

#Need Newton polynomial for newton model
#Can we find the elliptic curve of Lattes?

my_cursor.execute("""
CREATE TABLE families_dim_1 (
    label varchar(%s) PRIMARY KEY,
    degree integer,
    num_parameters integer,
    base_field base_field_type,
    sigma_invariants sigma_invariants_type,
    citations integer[],
    original_model model_type,
    is_polynomial boolean,
    is_Newton boolean,
    is_Lattes boolean,
    num_critical_points integer,
    automorphism_group_cardinality integer
    )
""",[function_label_length])


my_session.commit()
