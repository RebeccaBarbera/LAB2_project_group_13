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
the positive dataset was retrieved with the following query:
*positive_url = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment:false%29 AND (taxonomy_id:2759) AND (length:[40 TO ]) AND (reviewed:true) AND (existence:1) AND (ft_signal_exp:)%29&size=500"*
##### Negative dataset:
- absence of SP sequence
- experimental evidence for non SP-related compartments
##### Both positive and negative datasets:
- protein length
- protein evidence
- protein annotation status
- protein superkingdom
- Fragments

