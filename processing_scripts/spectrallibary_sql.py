import pandas as pd
from sqlite3 import connect

def get_metabolites(path_to_pos, path_to_neg):
    '''Create SQLite connection based on path to .db files''' 
    # make connection to local pos and neg ion-mode reference dbs
    pos_connection = connect(path_to_pos)
    neg_connection = connect(path_to_neg)

    # perform SQL queries to get dataframes of metabolites
    query = 'SELECT * FROM compondPrimary'
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
    # get unique inchikeys by ranking by ID
    unique = valid_inchikey[['ID', 'InChIKey']]
    unique.insert(loc=0, column='rank', 
                  value=unique.groupby(['InChIKey'])['ID'].rank())
    unique = unique.loc[unique['rank']==1]
    merged = unique.merge(chebi_to_inchikey, left_on='InChIKey', right_on='InChIKey')
    return merged[['InChIKey', 'chebi_id']]

def wraps_the_script(path_to_pos, path_to_neg, path_to_ref):
    chebi_to_inchikey = pd.read_csv(path_to_ref)
    annots = get_metabolites(path_to_pos, path_to_neg)
    return merge_dbs_for_conversion(annots, chebi_to_inchikey)