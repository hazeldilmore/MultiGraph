def maxquant_inference(path, experiment_id):
    '''Loads the maxquant proteinGroups file at the given file path and creates a relationship
    mapping protein ids to peptide ids. This is the raw mapping given by MaxQuant (is not well 
    validated or curated).'''
    df = pd.read_table(path)
    # ENH: check that df is of expected format
    
    # expand the protein and peptide ID mapping
    protein_ids = df['Protein IDs'].str.split(';', expand=True).stack().reset_index()
    peptide_seq = df['Peptide sequences'].str.split(';', expand=True).stack().reset_index()
    
    # rename columns to more accurately reflect what they are 
    protein_ids.rename(columns={0:'uniprot_id'}, inplace=True)
    peptide_seq.rename(columns={0:'peptide_sequence'}, inplace=True)
    
    # inner join on level_0 which represents original index of df
    merged = peptide_seq[['level_0', 'peptide_sequence']].merge(protein_ids[['level_0', 'uniprot_id']],
								right_on='level_0', left_on='level_0')
    merged['experiment_id'] = experiment_id
    return merged[['uniprot_id', 'peptide_sequence', 'experiment_id']]
