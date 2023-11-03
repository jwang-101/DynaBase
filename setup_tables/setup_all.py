"""
This file sets up the PostgreSQL tables for the
Database of Arithmetic Dynamics.
There are two tables for dimension 1 functions; one for number fields
and one for finite fields.
    functions_dim_1_NF
    functions_dim_1_FF

Auxilary Tables
    graphs_dim_1_NF
    rational_preperiodic_dim_1_NF

There are two tables for fields data
    number_fields
    finite_fields
    
There is one table for families of maps.

There is one table for citations.


Note that the function setup file contains custom types that
are used in subsequent tables.

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


# assumes that the path to root of the project has been added
# sage: load_attach_path('[path]')


###################################
###connect to database

load("connect.py")

#returns
# my_session - database connection
# my_cursor - cursor to my_session

#################################
###global constants

field_label_length = 15
function_label_length= 25

#################################
####



load("setup_tables/functions_dim_1_setup.py")
load("setup_tables/fields_setup.py")
load("setup_tables/families_setup.py")
#load("setup_tables/citations_setup.py")

my_session.commit()

