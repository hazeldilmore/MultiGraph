// Create nodes of the graph

// metabolite nodes
LOAD CSV WITH HEADERS FROM 'file:///metabolites.csv' AS row
MERGE (metabolite:Metabolite {metaboliteID: row.metabolite_id})
ON CREATE SET metabolite.name = row.metabolite_name, metabolite.chemical_formula = row.chemical_formula, 
metabolite.molecular_weight = toFloat(row.molecular_weight); 

// reaction nodes
LOAD CSV WITH HEADERS FROM 'file:///reactions.csv' AS row
MERGE (reaction:Reaction {reactionID: row.reaction_id})
ON CREATE SET reaction.name = row.reaction_name;

// pathway nodes
LOAD CSV WITH HEADERS FROM 'file:///pathways.csv' AS row
MERGE (pathway:Pathway {pathwayID: row.pathway_id})
ON CREATE SET pathway.name = row.pathway_name;

// protein nodes
LOAD CSV WITH HEADERS FROM 'file:///proteins.csv' AS row
MERGE (protein:Protein {proteinID: row.protein_id})
ON CREATE SET protein.product = row.product_name, protein.gene = row.gene_name;

// Create indexes on the data 
CREATE INDEX metabolite_id FOR (m:Metabolite) ON (m.metaboliteID);
CREATE INDEX metabolite_name FOR (m:Metabolite) ON (m.name);
CREATE INDEX reaction_id FOR (r:Reaction) ON (r.reactionID);
CREATE INDEX reaction_name FOR (r:Reaction) ON (r.name);
CREATE INDEX pathway_id FOR (p:Pathway) ON (p.pathwayID);
CREATE INDEX pathway_name FOR (p:Pathway) ON (p.name);
CREATE INDEX protein_id FOR (p:Protein) ON (p.proteinID);
CREATE INDEX protein_product FOR (p:Protein) ON (p.product);
CREATE INDEX protein_gene FOR (p:Protein) ON (p.gene);
CALL db.awaitIndexes();

// Creating relationships between the nodes 

// reaction to metabolites
LOAD CSV WITH HEADERS FROM 'file:///reaction_to_metabolite.csv' AS row
MATCH (reaction:Reaction {reactionID: row.reaction_id})
MATCH (metabolite:Metabolite {metaboliteID: row.metabolite_id})
MERGE (reaction)-[rm:CONTAINS]->(metabolite)
    ON CREATE SET rm.role = row.role;

// pathway to reaction 
LOAD CSV WITH HEADERS FROM 'file:///pathway_to_reaction.csv' AS row
MATCH (pathway:Pathway {pathwayID: row.pathway_id})
MATCH (reaction:Reaction {reactionID: row.reaction_id})
MERGE (pathway)-[pr:INCLUDES]->(reaction);

// something wrong with this one because it doesn't create any relationships in the graph
// i think the problem was that there was a NaN in the index column 
LOAD CSV WITH HEADERS FROM 'file:///protein_to_reaction.csv' AS row
MATCH (reaction:Reaction {reactionID: row.reaction_id})
MATCH (protein:Protein {proteinID: row.protein_id})
MERGE (protein)-[rp:CATALYZES]->(reaction);