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

def maxquant_mapping(path):
    '''Loads the maxquant proteinGroups file at the given file path and creates a relationship
    mapping sample names to uniprot IDs, and keep LFQ intensity as an attribute. This is not the 
    ideal relationship (in the future, we will have another inference engine translate this, not 
    MaxQuant).'''
    df = pd.read_table(path)
    # ENH check that df is of expected path 
    
    # expand the protein ID mapping
    protein_ids = df['Protein IDs'].str.split(';', expand=True).stack().reset_index()
    protein_ids.rename(columns={0:'uniprot_id'}, inplace=True)
    
    # pull out the unique sample names in this experiment 
    sample_names = [name[9:] for name in df.columns[df.columns.str.contains('Peptides ')]]
    
    # set up a dataframe 
    final_cols = ['uniprot_id', 'lfq_intensity', 'sample_name']
    out_df = pd.DataFrame(columns=final_cols)
    
    # loop through each sample file and add results to output dataframe 
    for name in sample_names: 
        # get peptide <-> protein mapping with nonzero LFQ intensity 
        nonzero = df.loc[df['LFQ intensity ' + name] != 0]
        
        # inner join between protein_ids and the nonzero LFQ intensities
        merged = protein_ids.merge(nonzero[['LFQ intensity ' + name]], left_on='level_0', right_index=True)
        # rename the column to a consistent 'LFQ intensity'
        merged.rename(columns={'LFQ intensity ' + name: 'lfq_intensity'}, inplace=True)
        # add sample name as a column                            
        merged['sample_name'] = name
        # concatenate to output dataframe 
        out_df = pd.concat([out_df, merged[final_cols]])
    return out_df
