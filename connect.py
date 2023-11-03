###################################
###connect to database

import psycopg2
import psycopg2.extras
try:
    my_session = psycopg2.connect(database="dad",
                        host="localhost",
                        user="dad_user",
                        password="dad_pass",
                        port="5432")
except:
    print("Unable to connect to the database") 
    
my_cursor = my_session.cursor(cursor_factory=psycopg2.extras.DictCursor)


