// Create nodes of the graph

// Metabolites nodes 
LOAD CSV WITH HEADERS FROM 'file:///metabolites.csv' AS row
MERGE (metabolite:Metabolite {ramp_id: row.ramp_id});
// is this correct when there aren't any other attributes? 

// Synonym nodes
LOAD CSV WITH HEADERS FROM 'file:///synonyms.csv' AS row
MERGE (synonym:Synonym {analyte_synonym_id: row.analyte_synonym_id})
ON CREATE SET synonym.synonym = row.synonym; 

// Chem_Props nodes
LOAD CSV WITH HEADERS FROM 'file:///chem_props.csv' AS row
MERGE (chem_props:Chem_Props {chem_source_id: row.chem_source_id})
ON CREATE SET chem_props.inchi = row.inchi, chem_props.inchi_key = row.inchi_key, 
chem_props.inchi_key_prefix = row.inchi_key_prefix, chem_props.iso_smiles = row.iso_smiles,
chem_props.mw = toFloat(row.mw), chem_props.monoisotop_mass = toFloat(row.monoisotop_mass), 
chem_props.common_name = row.common_name, chem_props.mol_formula = row.mol_formula;

// Class nodes
LOAD CSV WITH HEADERS FROM 'file:///metabolite_class.csv' AS row
MERGE (metabolite_class:Class {metabolite_class_id: row.meatbolite_class_id})
ON CREATE SET metabolite_class.class_source_id = row.class_source_id, 
metabolite_class.class_level_name = row.class_level_name, metabolite_class.class_name = row.class_name, 
metabolite_class.source = row.source; 

// Pathway nodes
LOAD CSV WITH HEADERS FROM 'file:///pathways.csv' AS row
MERGE (pathway:Pathway {pathwayRampId: row.pathwayRampId})
ON CREATE SET pathway.sourceId = row.sourceId, pathway.type = row.type, pathway.pathwayCategory = row.pathwayCategory, 
pathway.pathwayName = row.pathwayName; 

// Metabolite Ontology nodes
LOAD CSV WITH HEADERS FROM 'file:///ontology.csv' AS row
MERGE (ontology:Ontology {rampOntologyId: row.rampOntologyId})
ON CREATE SET ontology.commonName = row.commonName, ontology.HMDBOntologyType = row.HMDBOntologyType, 
ontology.metCount = row.metCount; 

// Reference Feature nodes
LOAD CSV WITH HEADERS FROM 'file:///features.csv' AS row
MERGE (reference_feature:Reference_Feature {reference_feature_id: row.reference_feature_id})
ON CREATE SET reference_feature.name = row.name, reference_feature.formula = row.formula, 
reference_feature.precursor_mz = row.precursor_mz, reference_feature.inchikey = row.inchikey, 
reference_feature.retention_time = row.retention_time; 

// Sample nodes
LOAD CSV WITH HEADERS FROM 'file:///samples.csv' AS row
MERGE (sample:Sample {sample_id:row.sample_id})
ON CREATE SET sample.filename = row.filename, sample.condition = row.condition, sample.replicate = row.replicate, 
sample.ion_mode = row.ion_mode; 

// Experiment nodes 
LOAD CSV WITH HEADERS FROM 'file:///experiments.csv' AS row
MERGE (experiment:Experiment {experiment_id:row.experiment_id})
ON CREATE SET experiment.experiment_desc = row.experiment_desc, experiment.mass_spec = row.mass_spec,
experiment.chromatography = row.chromatography; 

// Create indexes on the data 
CREATE INDEX ramp_id FOR (m:Metabolite) ON (m.ramp_id);
CREATE INDEX metabolite_synonym_id FOR (s:Synonym) ON (s.analyte_synonym_id);
CREATE INDEX chem_props_id FOR (c:Chem_Props) ON (c.chem_source_id);
CREATE INDEX metabolite_class_id FOR (c:Class) ON (c.meatbolite_class_id);
CREATE INDEX pathway_ramp_id FOR (p:Pathway) ON (p.pathwayRampId);
CREATE INDEX ontology_ramp_id FOR (o:Ontology) ON (o.rampOntologyId);
CREATE INDEX reference_feature_id FOR (r:Reference_Feature) ON (r.reference_feature_id);
CREATE INDEX sample_id FOR (s:Sample) ON (s.sample_id);
CREATE INDEX experiment_id FOR (e:Experiment) ON (e.experiment_id);
CALL db.awaitIndexes();

// Create relationships of the graph 
// metabolite to class 
LOAD CSV WITH HEADERS FROM 'file:///has_class.csv' AS row
MATCH (metabolite:Metabolite {ramp_id: row.ramp_id})
MATCH (metabolite_class:Class {metabolite_class_id: row.metabolite_class_id})
MERGE (metabolite)-[rm:HAS_CLASS]->(metabolite_class);

// pathway to metabolite 
LOAD CSV WITH HEADERS FROM 'file:///contains_metabolite.csv' AS row
MATCH (pathway:Pathway {pathwayRampId: row.pathwayRampId})
MATCH (metabolite:Metabolite {ramp_id: row.ramp_id})
MERGE (pathway)-[pr:CONTAINS_METABOLITE]->(metabolite)
    ON CREATE SET pr.pathwaySource = row.pathwaySource;

// metabolite to ontology
LOAD CSV WITH HEADERS FROM 'file:///has_ontology.csv' AS row
MATCH (metabolite:Metabolite {ramp_id: row.ramp_id})
MATCH (ontology:Ontology {rampOntologyId: row.rampOntologyId})
MERGE (metabolite)-[mo:HAS_ONTOLOGY]->(ontology); 

// metabolite to synonym
LOAD CSV WITH HEADERS FROM 'file:///has_synonym.csv' AS row
MATCH (metabolite:Metabolite {ramp_id: row.ramp_id})
MATCH (synonym:Synonym {analyte_synonym_id: row.analyte_synonym_id})
MERGE (metabolite)-[ms:HAS_SYNONYM]->(synonym)
    ON CREATE SET ms.source = row.source;

// metabolite to chem_props
LOAD CSV WITH HEADERS FROM 'file:///has_properties.csv' AS row
MATCH (metabolite:Metabolite {ramp_id: row.ramp_id})
MATCH (chem_props:Chem_Props {chem_source_id: row.chem_source_id})
MERGE (metabolite)-[mc:HAS_PROPERTIES]->(chem_props)
    ON CREATE SET mc.chem_data_source = row.chem_data_source;

// reference feature to metabolite
LOAD CSV WITH HEADERS FROM 'file:///matches_to.csv' AS row
MATCH (reference_feature:Reference_Feature {reference_feature_id: row.reference_feature_id})
MATCH (metabolite:Metabolite {ramp_id: row.ramp_id})
MERGE (reference_feature)-[rm:MATCHES_TO]->(metabolite);

// sample to experiment 
LOAD CSV WITH HEADERS FROM 'file:///part_of.csv' AS row
MATCH (sample:Sample {sample_id: row.sample_id})
MATCH (experiment:Experiment {experiment_id: row.experiment_id})
MERGE (sample)-[se:MATCHES_TO]->(experiment);