"""
Setup the PostgreSQL tables to work with number fields
and finite fields.

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
    DROP TABLE IF EXISTS finite_fields
""")


############################
# create custom types


#create new types


######################################
# Create Fields table Schema

#do we need embeddings?             F['embeddings']
# label = NF.label, FF.label, pF.label ?

my_cursor.execute("""
CREATE TABLE finite_fields (
    label varchar(%s) PRIMARY KEY,
    modulus_coeffs integer[],
    characteristic integer,
    cardinality integer
  )
""",[field_label_length])

my_session.commit()
