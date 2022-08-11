import pandas as pd

def read_and_filter_msdial_file(path_to_alignment):
    '''Read .msdial file that is output from MS-DIAL pipeline on katherine-johnson and filter
    out metabolites that do not have an INCHIKEY or have an identification that is not based
    on MS2 spectra.'''
    df = pd.read_csv(path_to_alignment, sep='\t', skiprows=4, index_col='Alignment ID')
    return df.loc[(df['INCHIKEY'].notna()) & (~df['Metabolite name'].str.contains('w/o MS2'))]

def get_combined_annots(path_to_positive, path_to_negative, path_to_chebi):
    '''Combine the positive and negative annotations'''
    cols_of_interest = ['Average Rt(min)', 'Average Mz', 'Metabolite name', 'Adduct type', 'Post curation result',
                         'Formula', 'INCHIKEY']
    positive_aligned = read_and_filter_msdial_file(path_to_positive)[cols_of_interest]
    negative_aligned = read_and_filter_msdial_file(path_to_negative)[cols_of_interest]
    combined_aligned = pd.concat([positive_aligned, negative_aligned])
    
    # read file with InChIKey to ChEBI mapping 
    chebi_to_inchikey = pd.read_csv(path_to_chebi)
    
    return combined_aligned.merge(chebi_to_inchikey, left_on='INCHIKEY', right_on='InChIKey', how='left')

def wraps_the_script(path_to_positive, path_to_negative, path_to_chebi):
    combined = get_combined_annots(path_to_positive, path_to_negative, path_to_chebi) 
    
    # print some warning message about all the ones where chebi_id.isna() 
    
    # for those where chebi_id.notna(), upload to graph database with relationship info