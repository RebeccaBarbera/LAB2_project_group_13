# LAB2_project_group_13

## Signal peptide prediction
The aim of this project is to evaluate and compare different computational methods for detecting signal peptides as well as addressing the subproblem of subcellular localisation and protein function prediction. 

## Software and tools needed
- `Python 3` → main programming language for data processing.
- `Biopython (Bio.SeqIO)` → for handling FASTA input/output.
- `Requests` → for making HTTP requests to UniProt REST API.
- `GitHub` / `Git` → for version control and collaboration

### 1. Data collection
The first step is to retreive both positive and negative dataset for evaluation.
##### database used: `https://www.uniprot.org`

#### Common criteria for protein selection:
##### Both positive and negative datasets:
- protein length
- protein evidence
- protein annotation status
- protein superkingdom
- Fragments
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
To ensure these proteins lacked signal peptides, proteins were chosen from experimentally validated subcellular localizations that are not routed through the secretory pathway, including cytoplasm (SL-0091), nucleus (SL-0191), mitochondrion (SL-0173), plastid (SL-0209), chloroplast (SL-0204), and cytoskeleton (SL-0039):
- negative_url = `"https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment:false%29 AND (reviewed:true) AND (existence:1) AND (length:[40 TO ]) AND (taxonomy_id:2759) NOT (ft_signal:) AND ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204) OR (cc_scl_term_exp:SL-0039))%29&size=500"`
##### Results from both positive and negative datasets where retrieved in JSON format. 

### 2. Data processing pipeline
#### the next step is to filter the positive dataset to meet the following criteria:
- No fragments
- Select only eukaryotic proteins
- Filter-out sequences shorter than 40 residues
- Filter-out unreviewed proteins
- Select on protein with experimental SP evidence
- Filter out proteins with SP shorter than 14 residues
- Protein existence: evidence at protein level
- Existence of the cleavage site
##### To filter the positive dataset the custom python scrypt named `positive_set.py` was used.
##### The output file is the `positive_filtered.tsv` file with the following information:
1. The protein UniProt accession number
2. The organism's name
3. The Eukaryotic kingdom (Metazoa, Fungi, Plants, Other)
4. The protein length
5. The position of the signal peptide cleavage site

