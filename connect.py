###################################
###connect to database

import psycopg2
import psycopg2.extras

from config import load_config

config = load_config()
try:
    # connecting to the PostgreSQL server
    my_session = psycopg2.connect(**config)
    print('Connected to the PostgreSQL server.')

except (psycopg2.DatabaseError, Exception) as error:
    print(error)


my_cursor = my_session.cursor(cursor_factory=psycopg2.extras.DictCursor)


