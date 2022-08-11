import pandas as pd
import pymysql
from sqlalchemy import create_engine

def make_connection_str(db, user, password):
    '''Makes connection str for given database.'''
    return 'mysql+pymysql://' + user + ':' + password + '@localhost/' + db

def sql_to_dataframe(query, connection):
    '''Runs SQL query on given db connection, outputs pd.DataFrame'''
    return pd.read_sql(query, connection)

def get_inchikeys(user, password):
    '''Extract InChIKey to ChEBI dataframe through chebi's structures
    relational database.'''
    # make connection to locally installed chebi relational database
    chebi_connection_str = make_connection_str('chebi', user, password)
    chebi_connection = create_engine(chebi_connection_str)

    # select the InChIKeys
    query = "select * from structures where type='InChIKey'"
    chebi_to_inchikey = sql_to_dataframe(query, chebi_connection)
    chebi_to_inchikey.rename(columns={'compound_id':'chebi_id', 'structure':'InChIKey'}, 
                             inplace=True)
    return chebi_to_inchikey[['chebi_id', 'InChIKey']]

def chebi_to_inchikey(user, password):
    get_inchikeys(user, password).to_csv('chebi_inchikey_conversion.csv', index=False)