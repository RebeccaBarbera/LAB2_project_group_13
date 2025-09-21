# LAB2_project_group_13
## Signal peptide prediction
The aim of this project is to evaluate and compare different computational methods for detecting signal peptides as well as addressing the subproblem of subcellular localisation and protein function prediction. 

### 1. Data cllection
The first step is to retreive both positive and negative dataset for evaluation.

#### Common criteria for protein selection:
##### Positive dataset:
- signal peptide evidence 
- knowledge of Signal Peptide (SP) cleavage site
- SP length > 13
##### The positive dataset was retrieved from UniProt using a query that selected non-fragment, reviewed proteins from Eukaryota (taxonomy ID 2759) with a sequence length of at least 40 amino acids, evidence at protein level, and an experimentally annotated signal peptide
- positive_url = `"https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment:false%29 AND (taxonomy_id:2759) AND (length:[40 TO ]) AND (reviewed:true) AND (existence:1) AND (ft_signal_exp:)%29&size=500"`
##### Negative dataset:
- absence of SP sequence
- experimental evidence for non SP-related compartments
##### The negative dataset was retrieved from UniProt by selecting non-fragment, reviewed proteins from Eukaryota (taxonomy ID 2759) with evidence at the protein level and a minimum sequence length of 40 amino acids. 
To ensure these proteins lacked signal peptides, entries annotated with signal sequences were excluded. Instead, proteins were chosen from experimentally validated subcellular localizations that are not routed through the secretory pathway, including cytoplasm (SL-0091), nucleus (SL-0191), mitochondrion (SL-0173), plastid (SL-0209), chloroplast (SL-0204), and cytoskeleton (SL-0039):
- negative_url = `"https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment:false%29 AND (reviewed:true) AND (existence:1) AND (length:[40 TO ]) AND (taxonomy_id:2759) NOT (ft_signal:) AND ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039))%29&size=500"`
##### Both positive and negative datasets:
- protein length
- protein evidence
- protein annotation status
- protein superkingdom
- Fragments

