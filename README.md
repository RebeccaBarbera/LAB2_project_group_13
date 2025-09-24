# LAB2_project_group_13
#### Rebecca Barbera, Aniello Di Vaio, Nina Talajic, Domenico Zianni

## Signal peptide prediction
The aim of this project is to evaluate and compare different computational methods for detecting signal peptides as well as addressing the subproblem of subcellular localisation and protein function prediction. 
## Table of Contents
- [Software and tools needed](#software-and-tools-needed)
- [Data collection](#1-data-collection)
- [Data filtering pipeline](#2-data-filtering-pipeline)
  - [Output files for data collection](#output-files)
  - [Dataset summary table](#dataset-summary)
- [Data pre-processing](#data-pre-processing)
   - [Clustering](#clustering)
   - [Filtering into a TSV file](#filtering_into_a_TSV_file)
   - [Data summary table](#data_summary_table)
   - [Data split](#data_split)
     - [Data split overall results](#Data_split_overall_results)
  - [Five-fold Cross validation](#five_fold_cross_validation)

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
**Please note:** that the negative dataset was directly retrieved from UniProt using the query criteria, without the need for further filtering of the JSON response. The script is used only to extract the required fields and to format the results into TSV and FASTA files.

##### Both `positive_set.fasta` and `negative_set.fasta` are in standard FASTA format, where each entry begins with '>' followed by the UniProt accession and the following line contains the full amino acid sequence.

## Dataset Summary

| Dataset  | Total | Metazoa | Fungi | Viridiplantae | Other | N-terminal TM helix |
|----------|-------|---------|-------|---------------|-------|----------------------|
| Positive |  2932 |    2420 |   165 |           311 |    36 | -                    |
| Negative | 20615 |   12419 |  3727 |          4111 |   358 | 2477                 |


## 3. Data pre-processing
- The first step of data pre-processing consists in using clustering methods to remove non-reduntant sequences from the dataset.
- Next clustered data will be further split into teo sets:
  - the **`training set`**: used to train the methods, optimize model hyperparameters and perform
cross-validation experiments
  - the **`benchmark set`** (also known as the holdout set):  used to test the generalization performance of the different models
### Clustering
Clustering is executed with a software suite called **`MMseq2`**.
MMSeq2 is the fastest method available for clustering, due to its implementation of three distinct clustering modes:_ Greedy set cover, Greedy incremental, and Connected-component clustering._

**The following commands have been used to cluster both positive and negative datasets into 2 different clustered sets:**

For the positive datase:
- `mmseqs easy-cluster positive_set.fasta pos_cluster-results tmp_pos --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1`

For the negative dataset 
- `mmseqs easy-cluster negative_set.fasta neg_cluster-results tmp_neg --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1`

These commands take all sequences in *_set.fasta, compare them to each other and group them into clusters of similar sequences with ≥30% identity and ≥40% coverage. They save the results in *_cluster-results, and use _tmp_*_ as a temporary working directory for the program.

**Please note**: Prior to clustering, we converted the FASTA files from DOS to Unix format using **`dos2unix`**. This step was necessary because the original files contained trailing spaces at the end of sequences due to Windows formatting.

**Output files:**
**Positive dataset:**(N=x)
- `pos_cluster-results_all_seqs.fasta`
- `pos_cluster-results_cluster.tsv`
- **`pos_cluster-results_rep_seq.fasta`**

**Negative dataset**(N=x):
- `neg_cluster-results_all_seqs.fasta`
- `neg_cluster-results_cluster.tsv`
- **`neg_cluster-results_rep_seq.fasta`**

### Filtering into a TSV file
Both `pos_cluster-results_rep_seq.fasta` and `neg_cluster-results_rep_seq.fasta` were used to retrieve **representative sequences from both clusters** and extract sequence infomration (Kingdom, protein length, etc) from the original tsv file obtained from both the original `positive_set.tsv and `negative_set.tsv`. 

To obtain the desired results **bash shell scripting was used accordingly:**

For the positive dataset:
- `grep "^>" pos_cluster-results_rep_seq.fasta | sed 's/^>//; s/[[:space:]]*$//' > positive_ids.txt`
- `head -n 1 positive_set.tsv > positive_info.tsv`
- `grep -F -f positive_ids.txt positive_set.tsv >> positive_info.tsv`

For the negative dataset:
- `grep "^>" neg_cluster-results_rep_seq.fasta | sed 's/^>//; s/[[:space:]]*$//' > negative_ids.txt`
- `head -n 1 negative_set.tsv > neg_info.tsv`
- `grep -F -f pnegative_ids.txt positive_set.tsv >> neg_info.tsv`

### Data summary table
| Dataset  | Total | Metazoa | Fungi | Viridiplantae | Other | N-terminal TM helix |
|----------|-------|---------|-------|---------------|-------|----------------------|
| Positive |  1092 |    866 |   95 |           103 |    28 | -                    |
| Negative | 8934 |   4697 |  2475 |          1594 |   168 | 900                 |

### Data split
The next step is to split the data into a 80/20 ratio, where **80%** belongs to the **training set** and the remaining **20%** belongs to the **benchmarking set**. This step is crucial for ensuring unbiased results and that the model learns generalizable patterns. 

**Extract IDs from the representative Fasta files obtained from teh MMSeq run**
- `grep "^>" pos_cluster-results_rep_seq.fasta | sed 's/^>//' > pos_ids.txt`
- `grep "^>" neg_cluster-results_rep_seq.fasta | sed 's/^>//' > neg_ids.txt`

**Shuffle IDs** 
- `sort -R pos_ids.txt > pos_shuffled_ids.txt`
- `sort -R neg_ids.txt > neg_shuffled_ids.txt`

**Calculate 80% split sizes**
- `pos_total=$(wc -l < pos_shuffled_ids.txt)`
- `pos_train_lines=$(( pos_total * 80 / 100 ))`

- `neg_total=$(wc -l < neg_shuffled_ids.txt)`
- `neg_train_lines=$(( neg_total * 80 / 100 ))`

**Split into training / benchmarking ID lists**
- `head -n $pos_train_lines pos_shuffled_ids.txt > pos_train_ids.txt`
- `tail -n +$((pos_train_lines+1)) pos_shuffled_ids.txt > pos_benchmark_ids.txt`

- `head -n $neg_train_lines neg_shuffled_ids.txt > neg_train_ids.txt`
- `tail -n +$((neg_train_lines+1)) neg_shuffled_ids.txt > neg_benchmark_ids.txt`

**Extract FASTA sequences with Python script**
- `python3 get_seq.py pos_train_ids.txt pos_cluster-results_rep_seq.fasta pos_train.fasta`
- `python3 get_seq.py pos_benchmark_ids.txt pos_cluster-results_rep_seq.fasta pos_benchmark.fasta`

- `python3 get_seq.py neg_train_ids.txt neg_cluster-results_rep_seq.fasta neg_train.fasta`
- `python3 get_seq.py neg_benchmark_ids.txt neg_cluster-results_rep_seq.fasta neg_benchmark.fasta`

**Merge positives + negatives**
- `cat pos_train.fasta neg_train.fasta > train.fasta`
- `cat pos_benchmark.fasta neg_benchmark.fasta > benchmark.fasta`

#### Data split overall results
| Set       | Positive | Negative | Total |
|-----------|----------|----------|-------|
| Training  | 873      | 7147     | 8020  |
| Benchmark | 219      | 1787     | 2006  |

### Five-fold Cross Validation 
This step is to randomly split the training set into 5 different subsets, preserving the overall positive/negative ratio on each subset.

**First, we extract the IDs from both the merged training.fasta dataset and benchmark.fasta dataset:**
- grep "^>" train.fasta| sed 's/^>//' > train_ids.txt
- grep "^>" benchmark.fasta| sed 's/^>//' > bench_ids.txt
  
**Then we randomly shuffle the IDs**
- sort -R train_ids.txt > train_ids_shuffled.txt
- sort -R bench_ids.txt > bench_ids_shuffled.txt

**Lastly we split the dataset into 5 roughly equal folds**
- gsplit -n l/5 train_ids_shuffled.txt fold_
- gsplit -n l/5 bench_ids_shuffled.txt fold_bench_

#### output files:

**Training set:**
- `fold_aa`
- `fold_ab`
- `fold_ac`
- `fold_ad`
- `fold_ae`

**Benchmark set:**
- `fold_bench_aa`
- `fold_bench_ab`
- `fold_bench_ac`
- `fold_bench_ad`
- `fold_bench_ae`

