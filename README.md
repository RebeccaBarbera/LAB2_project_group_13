# LAB2_project_group_13
#### Rebecca Barbera, Aniello Di Vaio, Nina Talajic, Domenico Zianni

## Signal peptide prediction
The aim of this project is to evaluate and compare different computational methods for detecting signal peptides as well as addressing the subproblem of subcellular localisation and protein function prediction. 
## Table of Contents
- [Software and tools needed](#software-and-tools-needed)
- [Data collection](#1-data-collection)
- [Data filtering pipeline](#2-data-filtering-pipeline)
  - [Output files for data collection](#output-files)
- [Data pre-processing](#data-pre-processing)
   - [Clustering](#clustering)

## Software, pakcages and tools needed
- `Python 3` → main programming language for data processing.
- `Biopython (Bio.SeqIO)` → for handling FASTA input/output.
- `Requests` → for making HTTP requests to UniProt REST API.
- `GitHub` / `Git` → for version control and collaboration.
- `MMSeqs2` → software suite used for clustering.

## 1. Data collection
The first step is to retrieve both positive and negative dataset for evaluation.
- database used: `https://www.uniprot.org`

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

## 2. Data filtering pipeline
#### The next step is to filter the dataset to meet the following criteria:
##### Positive dataset: 
- No fragments
- Select only eukaryotic proteins
- Filter-out sequences shorter than 40 residues
- Filter-out unreviewed proteins
- Select only proteins with experimental SP evidence
- Filter out proteins with SP shorter than 14 residues
- Protein existence: evidence at protein level
- Existence of the cleavage site
##### Negative dataset:
- No fragments
- Filter-out unreviewed proteins
- Protein existence: evidence at protein level
- Select only eukaryotic proteins
- Filter-out sequences shorter than 40 residues
- Filter-out sequences having SP (any evidence)
- Select only proteins experimentally verified to be localized into: cytosol, nucleus, mitochondrion, plastid, peroxisome, cell membrane.

### To filter both dataset the custom python script named `data-gathering.py` was used.[**](#neg-note)
##### Output files:
- `positive_set.tsv`
- `positive_set.fasta`
- `negative_set.tsv`
- `negative_set.fasta`
##### The `positive_set.tsv` and `negative_set.tsv` files contain the following information:
**For the `positive_set.tsv`:**
1. The protein UniProt accession number
2. The organism's name
3. The Eukaryotic kingdom (Metazoa, Fungi, Plants, Other)
4. The protein length
5. The position of the signal peptide cleavage site

**For the `negative_set.tsv`:**
1. The protein UniProt accession
2. The organism name
3. The Eukaryotic kingdom (Metazoa, Fungi, Plants, Other)
4. The protein length
5. Whether the protein has a transmembrane helix starting in the first 90 residues
(true or false)

<a name="neg-note"></a>
###### **Please note that the negative dataset was directly retrieved from UniProt using the query criteria, without the need for further filtering of the JSON response. The script is used only to extract the required fields and to format the results into TSV and FASTA files.

##### Both `positive_set.fasta` and `negative_set.fasta` are in standard FASTA format, where each entry begins with '>' followed by the UniProt accession and the following line contains the full amino acid sequence.

## 3. Data pre-processing
- The first step of data pre-processing consists in using clustering methods to remove non-reduntant sequences from the dataset.
- Next clustered data will be further split into teo sets:
  - the **`training set`**: used to train the methods, optimize model hyperparameters and perform
cross-validation experiments
  - the **`benchmark set`** (also known as the holdout set):  used to test the generalization performance of the different models
### Clustering
Clustering is executed with a software suite called `MMseq2`
