import pandas as pd
import pymysql
from sqlite3 import connect
from sqlalchemy import create_engine

def make_connection_str(db, user='root', password=''):
    '''Makes connection str for given database.'''
    return 'mysql+pymysql://' + user + ':' + password + '@localhost/' + db

def sql_to_dataframe(query, connection):
    '''Runs SQL query on given db connection, outputs pd.DataFrame'''
    return pd.read_sql(query, connection)

def get_inchikeys(user, password):
    '''Extract InChIKey to ChEBI dataframe through chebi's structures
    relational database.'''
    # make connection to locally installed chebi relational database
    chebi_connection_str = make_connection_str(db='chebi')
    chebi_connection = create_engine(chebi_connection_str)

    # select the InChIKeys
    query = "select * from structures where type='InChIKey'"
    chebi_to_inchikey = sql_to_dataframe(query, chebi_connection)
    chebi_to_inchikey.rename(columns={'compound_id':'chebi_id'}, inplace=True)
    return chebi_to_inchikey[['chebi_id', 'InChIKey']]

def get_metabolites(path_to_pos, path_to_neg):
    '''Create SQLite connection based on path to .db files''' 
    # make connection to local pos and neg ion-mode reference dbs
    pos_connection = connect(path_to_pos)
    neg_connection = connect(path_to_neg)

    # perform SQL queries to get dataframes of metabolites
    query = 'SELECT * FROM compoundPrimary'
    pos_df = sql_to_dataframe(query, pos_connection)
    pos_df['ion_mode'] = 'Positive'
    neg_df = sql_to_dataframe(query, neg_connection)
    neg_df['ion_mode'] = 'Negative'
    return pd.concat([pos_df, neg_df])

def merge_dbs_for_conversion(annots, chebi_to_inchikey):
    '''Concatenates reference library of metabolite annotations with 
    ChEBI to InChIKey reference to make a conversion between reference
    annotations and ChEBI IDs (therefore the Reactome db).'''
    valid_inchikey = annots.loc[(annots['InChIKey']!='') & 
                                (annots['InChIKey'] != 'Internal Standard')]
    return valid_inchikey.merge(chebi_to_inchikey, left_on='InChIKey', right_on='structure')

def make_mb_nodes_relationships(path_to_neg, path_to_pos, user='root', password=''):
    '''This function creates a .csv file corresponding to the MATCHES_TO relationship
    and a .csv corresponding to the :Metabolite node in the graph database.'''
    chebi_to_inchikey = get_inchikeys(user, password)
    all_metabolites = get_metabolites(path_to_pos, path_to_neg)
    merged_df = merge_dbs_for_conversion(all_metabolites, chebi_to_inchikey)
    
    # make .csv corresponding to :Metabolite node
    all_metabolites.to_csv('metabolites.csv')
    # make .csv corresponding to MATCHES_TO relationship
    merged_df.to_csv('MATCHES_TO.csv')