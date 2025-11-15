# Comparative Analysis of Eukaryotic Signal Peptide Predictors: Integrating the von Heijne PSWM and Support Vector Machines (SVM)

This repository gathers all resources and analyses developed within the Laboratory of Bioinformatics 2 course, part of the Master’s Programme in Bioinformatics at the University of Bologna.

This project was executed and completed by: [Rebecca Barbera](https://github.com/RebeccaBarbera), [Aniello Di Vaio](https://github.com/aniello-divaio), [Nina Talajic](https://github.com/nicnaaa), and [Domenico Zianni](https://github.com/mingo002).

---

## Signal peptide prediction
The aim of this project is to evaluate and compare different computational methods for detecting signal peptides as well as addressing the subproblem of subcellular localisation and protein function prediction. 

## Table of Contents

| | |
| :--- | :--- |
| **1.** [Software and tools needed](#software-and-tools-needed) | **4.** [Data Splitting](#data-splitting) |
| **2.** [Data Collection](#data-collection) | &nbsp;&nbsp;&nbsp; `4.1` [Data split overall results](#data-split-overall-results) |
| &nbsp;&nbsp;&nbsp; `2.1` [Output files for data collection](#output-files-for-data-collection) | &nbsp;&nbsp;&nbsp; `4.2` [Five-fold Cross Validation](#five-fold-cross-validation) |
| &nbsp;&nbsp;&nbsp; `2.2` [Dataset summary table](#dataset-summary-table) | **5.** [Data analysis](#data-analysis) |
| **3.** [Data Filtering Pipeline](#data-filtering-pipeline) | &nbsp;&nbsp;&nbsp; `5.1` [Comparative Amino Acid Composition](#comparative-amino-acid-composition) |
| &nbsp;&nbsp;&nbsp; `3.1` [Clustering](#clustering) | &nbsp;&nbsp;&nbsp; `5.2` [Taxonomic Classification](#taxonomic-classification) |
| &nbsp;&nbsp;&nbsp; `3.2` [Filtering into a TSV file](#filtering-into-a-tsv-file) | **6.** [Von Heijen Model](#von-heijen-model) |
| &nbsp;&nbsp;&nbsp; `3.3` [Data clustering table](#data-clustering-table) | **7.** [Support Vector Machine (SVM)](#support-vector-machine-svm) |
| | **8.** [Performance Evaluation](#performance-evaluation) |
## Software, packages and tools needed

> `Python 3` → main programming language for data processing.
> 
> `Biopython (Bio.SeqIO)` → for handling FASTA input/output.
>
> `Requests` → for making HTTP requests to UniProt REST API.
>
> `GitHub` / `Git` → for version control and collaboration.
>
> `MMSeqs2` → software suite used for clustering sequences.
>
> `NumPy` → for efficient numerical operations and handling multi-dimensional arrays (e.g., PSWM calculations).
>
> `Pandas` → for data manipulation, DataFrame management, and analyzing tabular data (TSV files).
>
> `Matplotlib` → static visualizations (density plots, heatmaps, confusion matrices).
>
> `Biopython (ProtParam)` → extraction of physicochemical features (`pI`, `MW`, `hydrophobicity`).
>
> `Scikit-learn` → ML pipeline:
> - `RandomForestClassifier` → feature ranking.
> - `SVC` → model training.
> - `StandardScaler` → feature normalization.
> - `metrics` → evaluating model performance.
>
> `Biopython` → Biological sequence analysis, FASTA parsing, and data handling for bioinformatics
>
> `Counter` → Counting and summarizing occurrences of amino acids or sequence elements
>
> `SeqIO` *(from Biopython)* → FASTA parsing and sequence extraction

## 1. Data collection
The first step is to retrieve both positive and negative dataset for evaluation. The database used for this purpose is [UniProt](https://www.uniprot.org).

#### Common criteria for protein selection:
- protein length
- protein evidence
- protein annotation status
- protein superkingdom
- Fragments
##### Positive dataset:
- signal peptide evidence 
- knowledge of Signal Peptide (SP) cleavage site
- SP length > 13

The positive dataset was retrieved from UniProt using a query that selected non-fragment, reviewed proteins from Eukaryota (taxonomy ID 2759) with a sequence length of at least 40 amino acids, evidence at protein level, and an experimentally annotated signal peptide. 

##### Negative dataset:
- absence of SP sequence
- experimental evidence for non SP-related compartments

The negative dataset was retrieved by selecting non-fragment, reviewed proteins from Eukaryota (taxonomy ID 2759) with evidence at the protein level and a minimum sequence length of 40 amino acids. 

To ensure these proteins lacked signal peptides, proteins were chosen from experimentally validated subcellular localizations that are not routed through the secretory pathway, including: cytoplasm (SL-0091), nucleus (SL-0191), mitochondrion (SL-0173), plastid (SL-0209), chloroplast (SL-0204) and cytoskeleton (SL-0039).

Results from both positive and negative datasets where retrieved in JSON format. 

### Data filtering pipeline

---

The next step is to filter the dataset to meet common criteria (Positive and Negative):
- No fragments
- Select only eukaryotic proteins
- Filter-out sequences shorter than 40 residues
- Filter-out unreviewed proteins
- Protein existence: evidence at protein level
and specific criteria per datase:
##### Positive Dataset Criteria
- Select only proteins with experimental SP evidence
- Filter out proteins with SP shorter than 14 residues
- Existence of the cleavage site
##### Negative Dataset Criteria
- Filter-out sequences having SP (any evidence)
- Select only proteins experimentally verified to be localized into: cytosol, nucleus, mitochondrion, plastid, peroxisome, cell membrane.

To filter both dataset the custom python script `data-gathering.py` was used.

`positive_set.tsv` and `negative_set.tsv` files was obtained, both containing these parameters:
- The protein UniProt accession number
- The organism's name
- The Eukaryotic kingdom (Metazoa, Fungi, Plants, Other)
- The protein length

A fifth parameter:
- For `positive_set.tsv`: The position of the signal peptide cleavage site.
- For `negative_set.tsv`: Whether the protein has a transmembrane helix starting in the first 90 residues (true or false).

<a name="neg-note"></a>
###### **Please note:** that the negative dataset was directly retrieved from UniProt using the query criteria, without the need for further filtering of the JSON response. The script is used only to extract the required fields and to format the results into TSV and FASTA files.

###### Both `positive_set.fasta` and `negative_set.fasta` are in standard FASTA format, where each entry begins with '>' followed by the UniProt accession and the following line contains the full amino acid sequence.

#### Dataset Summary

| Dataset  | Total | Metazoa | Fungi | Viridiplantae | Other | N-terminal TM helix |
|----------|-------|---------|-------|---------------|-------|----------------------|
| Positive |  2932 |    2420 |   165 |           311 |    36 | -                    |
| Negative | 20615 |   12419 |  3727 |          4111 |   358 | 2477                 |


## 2. Data pre-processing and Data Clustering
The first step of data pre-processing consists in using clustering methods to remove non-reduntant sequences from the dataset. Clustered data will be further split into two sets:
  - the **`training set`**: used to train the methods, optimize model hyperparameters and perform
cross-validation experiments
  - the **`benchmark set`** (also known as the holdout set):  used to test the generalization performance of the different models

**Clustering** is executed with a software suite called **`MMseq2`**, the fastest method available for clustering, due to the implementation of three distinct clustering modes: `Greedy set cover`, `Greedy incremental` and `Connected-component clustering`.

The following commands have been used to cluster both positive and negative datasets into two different clustered sets:

For positive:
- `mmseqs easy-cluster positive_set.fasta pos_cluster-results tmp_pos --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1`

For negative: 
- `mmseqs easy-cluster negative_set.fasta neg_cluster-results tmp_neg --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1`

These commands take all sequences in *_set.fasta, compare them to each other and group them into clusters of similar sequences with ≥30% identity and ≥40% coverage. They save the results in *_cluster-results, and use _tmp_*_ as a temporary working directory for the program.

###### **Please note**: Prior to clustering, we converted the FASTA files from DOS to Unix format using **`dos2unix`**. This step was necessary because the original files contained trailing spaces at the end of sequences due to Windows formatting.

#### Filtering into a TSV file
Both `pos_cluster-results_rep_seq.fasta` and `neg_cluster-results_rep_seq.fasta` were used to retrieve **representative sequences** from both clusters and extract sequence information (Kingdom, Protein length, etc) from the original tsv file obtained from both the original `positive_set.tsv` and `negative_set.tsv`. 

To obtain the desired results **bash shell scripting** was used accordingly:

For the positive dataset:
1. `grep "^>" pos_cluster-results_rep_seq.fasta | sed 's/^>//; s/[[:space:]]*$//' > positive_ids.txt`
2. `head -n 1 positive_set.tsv > positive_info.tsv`
3. `grep -F -f positive_ids.txt positive_set.tsv >> positive_info.tsv`

For the negative dataset:
1. `grep "^>" neg_cluster-results_rep_seq.fasta | sed 's/^>//; s/[[:space:]]*$//' > negative_ids.txt`
2. `head -n 1 negative_set.tsv > neg_info.tsv`
3. `grep -F -f pnegative_ids.txt positive_set.tsv >> neg_info.tsv`

#### Data clustering summary table
| Dataset  | Total | Metazoa | Fungi | Viridiplantae | Other | N-terminal TM helix |
|----------|-------|---------|-------|---------------|-------|----------------------|
| Positive |  1092 |    866 |   95 |           103 |    28 | -                    |
| Negative | 8934 |   4697 |  2475 |          1594 |   168 | 900                 |

## 3. Data spliting
The next step is to split the data into a 80/20 ratio, where **80%** belongs to the **training set** and the remaining **20%** belongs to the **benchmarking set**. This step is crucial for ensuring unbiased results and that the model learns generalizable patterns. 
```
# Extract IDs from the representative Fasta files obtained from teh MMSeq run
grep "^>" pos_cluster-results_rep_seq.fasta | sed 's/^>//' > pos_ids.txt
grep "^>" neg_cluster-results_rep_seq.fasta | sed 's/^>//' > neg_ids.txt

# Shuffle IDs
sort -R pos_ids.txt > pos_shuffled_ids.txt
sort -R neg_ids.txt > neg_shuffled_ids.txt

# Calculate 80% split sizes
## positive
pos_total=$(wc -l < pos_shuffled_ids.txt)
pos_train_lines=$(( pos_total * 80 / 100 ))

## negative
neg_total=$(wc -l < neg_shuffled_ids.txt)
neg_train_lines=$(( neg_total * 80 / 100 ))

# Split into training / benchmarking ID lists
head -n $pos_train_lines pos_shuffled_ids.txt > pos_train_ids.txt
tail -n +$((pos_train_lines+1)) pos_shuffled_ids.txt > pos_benchmark_ids.txt

head -n $neg_train_lines neg_shuffled_ids.txt > neg_train_ids.txt
tail -n +$((neg_train_lines+1)) neg_shuffled_ids.txt > neg_benchmark_ids.txt

# Extract FASTA sequences with Python script
python3 get_seq.py pos_train_ids.txt pos_cluster-results_rep_seq.fasta pos_train.fasta
python3 get_seq.py pos_benchmark_ids.txt pos_cluster-results_rep_seq.fasta pos_benchmark.fasta

python3 get_seq.py neg_train_ids.txt neg_cluster-results_rep_seq.fasta neg_train.fasta
python3 get_seq.py neg_benchmark_ids.txt neg_cluster-results_rep_seq.fasta neg_benchmark.fasta

# Merge positives + negatives
cat pos_train.fasta neg_train.fasta > train.fasta
cat pos_benchmark.fasta neg_benchmark.fasta > benchmark.fasta
```
#### Data split overall results
| Set       | Positive | Negative | Total |
|-----------|----------|----------|-------|
| Training  | 873      | 7147     | 8020  |
| Benchmark | 219      | 1787     | 2006  |

### Five-fold Cross Validation 
This step is to randomly split the training set into 5 different subsets, preserving the overall positive/negative ratio on each subset.
```
# Extract the IDs from both the merged training.fasta dataset and benchmark.fasta dataset:
grep "^>" train.fasta| sed 's/^>//' > train_ids.txt
grep "^>" benchmark.fasta| sed 's/^>//' > bench_ids.txt
  
# Randomly shuffle the IDs
sort -R train_ids.txt > train_ids_shuffled.txt
sort -R bench_ids.txt > bench_ids_shuffled.txt

# Split the dataset into 5 roughly equal folds
gsplit -n l/5 train_ids_shuffled.txt fold_
gsplit -n l/5 bench_ids_shuffled.txt fold_bench_
```

### Final TSV
The final TSV contained was organised in the following columns; the UniprotAccession code, Organism name ,Kingdom, Protein length, Signal Peptide Position, Positive/Negative set , and Fold Set.
The script used to generate the tsv file, `fold_tsv.ipynb` was written using pandas.

## 4. Data Analysis 
#### Distribution of Protein Lengths
Protein length distributions were visualized in R Studio using density plots for the positive and negative sequences in both the training and benchmark sets. To avoid distortion from a small number of very long sequences, the density plots were constructed on a logarithmic scale. 

---

#### Distrubution of SP position
The distribution of SP lengths were also visualized in R Studio using density plots for the positive sequences in both the training and benchmark sets. For each set, we plotted the distribution of SP lengths using density plots. The median SP length was calculated and displayed on the graph for both sets.

---

#### Comparative Amino Acid Composition
Amino acid composition of Signal Peptides (SPs) were compared against the background distribution of amino acids in SwissProt (data from [Expasy](https://web.expasy.org/docs/relnotes/relstat.html)). 
Extraction of all SP sequences, calculation of their amino acid frequencies, and plotting them against the SwissProt distribution were performed.

###### Taxonomic classification of the proteins was performed at both the kingdom and organism levels. The relative abundances of taxa in each dataset were visualized using pie charts.


#### Table of context
| Analysis                          | Visualization |
|-----------------------------------|---------------|
| Distribution of Protein Lengths    | [Protein Lengths distribution](4_data_analysis/Density_plot.png) |
| Distribution of SP Position        | [SPP distribution](4_data_analysis/SPPosition.png) |
| Comparative Amino Acid Composition | [AA Composition](4_data_analysis/AA_comparison.png) |
| Taxonomic Classification           | [Benchmark Classification](4_data_analysis/Kingdom_dist_bench.png) / [Training Classification](4_data_analysis/Kingdom_dist_train.png) |

## 5. The von Heijne Method
The main idea is to use use a ***Position-Speciﬁc Weight Matrix (PSWM)*** in order to model amino acid distribution around known cleavage sites. The retrieved scores were first stored in a ***PSPM*** and then a background model (Swiss-Prot database) was used as a reference amino acid distribution to compare our motifs against.

Detection of SP, comparison with Swiss-Prot and evaluation of the model were performed using the script [vonHejine.ipynb](5_vonHejine_model/vonHejine.ipynb)
 
 | *Fold* | *Training Folds* | *Validation Fold* | *Testing Fold* | *Motifs* | *Optimal Threshold* | *F1* | *Precision* | *Recall* | *MCC* |
| :------- | :----------------- | :------------------ | :--------------- | -----------: | --------------------: | -----: | ------------: | ---------: | ------: |
| 1        | ac, ad, ae         | ab                  | aa               |          531 |                 6.437 |  0.691 |         0.693 |      0.689 |   0.653 |
| 2        | aa, ad, ae         | ac                  | ab               |          528 |                 6.020 |  0.707 |         0.645 |      0.782 |   0.674 |
| 3        | aa, ab, ae         | ad                  | ac               |          519 |                 6.212 |  0.722 |         0.716 |      0.728 |   0.686 |
| 4        | aa, ab, ac         | ae                  | ad               |          522 |                 6.303 |  0.718 |         0.712 |      0.724 |   0.683 |
| 5        | ab, ac, ad         | aa                  | ae               |          519 |                 5.735 |  0.730 |         0.670 |      0.802 |   0.697 |
 
<table>
<tr>
<td>

<table>
<tr><th>Metric</th><th>Mean ± SE</th></tr>
<tr><td>F1 Score</td><td>0.714 ± 0.006</td></tr>
<tr><td>Precision</td><td>0.687 ± 0.012</td></tr>
<tr><td>Recall</td><td>0.745 ± 0.018</td></tr>
<tr><td>MCC</td><td>0.679 ± 0.007</td></tr>
<tr><td>Average Threshold</td><td>6.141 ± 0.121</td></tr>
</table>

</td>
<td>

<table>
<tr><th>Fold</th><th>TP</th><th>TN</th><th>FP</th><th>FN</th></tr>
<tr><td>fold_aa</td><td>122</td><td>1375</td><td>54</td><td>55</td></tr>
<tr><td>fold_ab</td><td>129</td><td>1370</td><td>71</td><td>36</td></tr>
<tr><td>fold_ac</td><td>131</td><td>1372</td><td>52</td><td>49</td></tr>
<tr><td>fold_ad</td><td>126</td><td>1380</td><td>51</td><td>48</td></tr>
<tr><td>fold_ae</td><td>142</td><td>1352</td><td>70</td><td>35</td></tr>
</table>

</td>
</tr>
</table>

[Heatmap](5_vonHejine_model/pswm_heatmap.png) and [Presicion Recall Curve](5_vonHejine_model/prc.png) were retrieved. 

## 6. Support Vector Machine
Building SVM models to enable a cross-model comparison between the Von Heijne method and machine learning–based prediction approaches. Extracting a comprehensive set of 31 physicochemical features from the N-terminal region of each protein sequence.
- **Amino Acid Composition** (20 features) computed over the first 22 residues (k=22).
- **Hydrophobicity** (3 features: max, mean, and std dev) based on the Kyte-Doolittle scale (k=40).
- **Secondary Structure** (3 features: helix, turn, and sheet fractions) predicted by ProteinAnalysis (k=40).
- **Charge Features** (2 features: Isoelectric Point and net charge at pH 7.0) (k=40).
- **Transmembrane Tendency** (2 features: max and mean) using the Zhao & London scale (k=40).
- **Refractivity** (2 features: max and mean) using the Jones D.D. scale (k=40).
  
The model is evaluated through 5-fold cross-validation using predefined splits.
In each fold, 3 subsets are used for training and 1 for validation. A Random Forest ranks the 31 features by Gini importance, and SVM accuracy is assessed using increasing top-k feature sets to select the optimal number of features. Grid Search (C, gamma) is performed on this optimal feature subset using the validation fold for scoring.

### Results per Fold ###
| Fold    |   Random Forest|   Accuracy |   Confusion Matrix |   Confusion Matrix Top Features |
|:--------|-----:|-----:|-----:|-----:|
| fold_1 |  [RF_Gini_fold1](6_SVM/RF_Gini_img/RF_Gini_fold1-2.png)| [Accuracy_fold1](6_SVM/ACC_plots/acc_vs_features_fold1.png) |   [CM_fold1](6_SVM/CMs_per_fold/CM_fold1_ALL.png) |   [CM_Top_fold1](6_SVM/CMs_per_fold/CM_fold_1_selected_feature.png) |
| fold_2 |  [RF_Gini_fold2](6_SVM/RF_Gini_img/RF_Gini_fold2-2.png) |  [Accuracy_fold2](6_SVM/ACC_plots/acc_vs_features_fold2.png) |   [CM_fold2](6_SVM/CMs_per_fold/CM_fold2_ALL.png) |   [CM_Top_fold2](6_SVM/CMs_per_fold/CM_fold_2_selected_feature.png) |
| fold_3 |  [RF_Gini_fold3](6_SVM/RF_Gini_img/RF_Gini_fold3-2.png) |  [Accuracy_fold3](6_SVM/ACC_plots/acc_vs_features_fold3.png) |   [CM_fold3](6_SVM/CMs_per_fold/CM_fold3_ALL.png) |   [CM_Top_fold3](6_SVM/CMs_per_fold/CM_fold_3_selected_feature.png) |
| fold_4 |  [RF_Gini_fold4](6_SVM/RF_Gini_img/RF_Gini_fold4-2.png) |  [Accuracy_fold4](6_SVM/ACC_plots/acc_vs_features_fold4.png) |   [CM_fold4](6_SVM/CMs_per_fold/CM_fold4_ALL.png) |   [CM_Top_fold4](6_SVM/CMs_per_fold/CM_fold_4_selected_feature.png) |
| fold_5 |  [RF_Gini_fold5](6_SVM/RF_Gini_img/RF_Gini_fold5-2.png) |  [Accuracy_fold5](6_SVM/ACC_plots/acc_vs_features_fold5.png) |   [CM_fold5](6_SVM/CMs_per_fold/CM_fold5_ALL.png) |   [CM_Top_fold5](6_SVM/CMs_per_fold/CM_fold_5_selected_feature.png) |

The MCC results for the cross-validation folds are shown in [MCC_per_fold](6_SVM/MCC_per_fold/MCC_PER_FOLD.png). 

## 7. Performance Evaluation
### The Vonhejine model

The VonHejine model performance was assessed using the [VH_bench.ipynb](7_performance_evalutation/VonHejine_perf/script/VH_bench.ipynb) script, which implement the full benchmarking pipeline. 

##### Methodology:

| **Step**                             | **Description**                                                                                                                                                                                                                              |
| ------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **1. Load Inputs**                   | Load the 20×15 PSWM (positions –13 to +2) and parse the benchmark FASTA into a dictionary of sequences.                                                                                                                                      |
| **2. Sliding-Window Scoring**        | Slide a 15-aa window over the first ~90 residues of each protein. For each window, sum the PSWM weights for all (aa, position) pairs. The highest-scoring window is stored as the predicted signal-peptide region and inferred –13 position. |
| **3. Ground-Truth Labels**           | Import metadata from `all.tsv` and map each protein to a binary label (SP vs non-SP).                                                                                                                                                        |
| **4. Threshold Classification**      | Apply a fixed cutoff (6.141) to convert model scores into predicted SP/non-SP classes. Generate TP, TN, FP, and FN counts.                                                                                                                   |
| **5. Performance Metrics**           | Compute F1, precision, recall, MCC, and visualize the confusion matrix.                                                                                                                                                                      |
| **6. Kingdom-Specific Analysis**     | Combine predictions with taxonomic metadata and compute TP/TN/FP/FN per kingdom. Visualize using full-outcome and error-only donut plots.                                                                                                    |
| **7. False-Negative Motif Analysis** | For false negatives, extract the true (–13 to +2) window relative to the experimentally annotated cleavage site and export as FASTA for sequence-logo analysis.                                                                              |

#### The confusion matrix summarizing the model’s predictions is as follows:

- **True Negatives (TN)**: 1548
- **False Positives (FP)**: 239
- **False Negatives (FN)**: 8
- **True Positives (TP**): 211

#### Based on this matrix, we obtain:

- **F1-score** = 0.631
- **Precision** = 0.469
- **Recall** = 0.963
- **MCC** = 0.620

### Kingdom-Level Percentage Distribution
To assess whether model performance varies across evolutionary groups, we computed per-kingdom TP, TN, FP, and FN statistics.
The following table reports the percentage contribution of each outcome within each kingdom:

| Kingdom        | TP (%) | TN (%) | FP (%) | FN (%) |
|----------------|--------|--------|--------|--------|
| Fungi          | 0.58   | 95.33  | 0.97   | 3.11   |
| Metazoa        | 3.00   | 81.70  | 2.73   | 12.57  |
| Other          | 1.53   | 83.16  | 2.55   | 12.76  |
| Viridiplantae  | 1.53   | 90.57  | 3.36   | 4.54   |

### Visualisation of outcome rates
#### Donut Plots by Kingdom:

| Kingdom        | Donut Plot |
|----------------|------------|
| Fungi          | [View Plot](7_performance_evalutation/VonHejine_perf/donuts/donutFungi.png) |
| Metazoa        | [View Plot](7_performance_evalutation/VonHejine_perf/donuts/donutMetazoa.png) |
| Other          | [View Plot](7_performance_evalutation/VonHejine_perf/donuts/donutOther.png) |
| Viridiplantae  | [View Plot](7_performance_evalutation/VonHejine_perf/donuts/donutViridiplantae.png) |

In the final stage of the evaluation, we built a [Sequence Logo](7_performance_evalutation/VonHejine_perf/logo.png) from the `false-negative` set to explore the motif characteristics that led the model to classify them incorrectly

### SVM
We next evaluated our SVM-based model using the following training and validation workflow in the [svm_bench.ipynb](7_performance_evalutation/SVM_perf/scripts/svm_bench.ipynb)

### SVM Benchmarking Workflow

| Step | Description |
|------|-------------|
| **1. Load trained SVM matrices** | Import the full training feature matrix (`training_features.npz`), including all engineered features, labels, and sequence IDs. Load the saved cross-validation results (`cv_results.npy`) containing best C, gamma, kernel, and top-k. |
| **2. Import benchmark sequences** | Parse the `benchmark.fasta` file into a dictionary of IDs → sequences. Load the true labels from `all.tsv` and map them to the benchmark IDs. |
| **3. Recreate feature extraction pipeline** | Recompute **identical features** used during training: (1) amino-acid composition, (2) hydrophobicity profiles, (3) secondary structure fractions, (4) pI and charge, (5) TM-tendency scale, (6) refractivity scale. Build a benchmark feature matrix with the same column order as the training matrix. |
| **4. Standardize benchmark features** | Apply the exact `StandardScaler` fitted on the full training set to ensure identical scaling for comparisons. |
| **5. Train final SVM models** | Train: (A) a baseline SVM using **all features**, and (B) the final SVM using only the **top-k features** determined by Random Forest Gini importance during training. |
| **6. Predict benchmark labels** | Apply both SVM models to the benchmark feature matrix to generate predictions on previously unseen sequences. |
| **7. Compute performance metrics** | Evaluate benchmark predictions using Accuracy, F1-score, MCC, and confusion matrices. Save confusion matrix images for both the all-feature and top-k models. |
| **8. Evaluate model sensitivity across taxa** | Map each benchmark sequence to its organism and quantify TP, TN, FP, and FN counts per organism. Generate stacked barplots to visualize species-level error distribution. |
| **9. Visualize global FP/FN distribution** | Produce whole-dataset pie charts for global FP/FN percentages. |
| **10. Species-level FP and FN exploration** | Generate side-by-side pie charts + barplots showing which organisms contribute most to FN and FP. |
| **11. Analyze false positives** | Test whether FP sequences disproportionately contain **N-terminal transmembrane helices** (via `TMHelixFirst90` metadata). Compute fraction and visualize which organisms contribute TM-driven false positives. |
| **12. Analyze false negatives** | Investigate why FN sequences fail: (i) compare N-terminal AA composition (train vs TP vs FN), (ii) compare SP lengths, (iii) compare top-k feature distributions. |
| **13. Export feature importance and diagnostic plots** | Save top-k selected features, MCC-vs-k curve, and all FN/FP species distributions for downstream interpretation. |

### Performance Comparison: All features vs. selected features SVM models

| Model                           | Accuracy | F1 Score | MCC   |
|---------------------------------|----------|----------|-------|
| **All Features Model**          | 0.974    | 0.884    | 0.870 |
| **Selected Features (Top-k)**   | 0.974    | 0.882    | 0.867 |

### Confusion matrix comparison: All features vs. selected features SVM models
| Model                           | TP | TN | FP   | FN|
|---------------------------------|----------|----------|-------|------|
| **All Features Model**          | 199    | 1755    | 32 | 20 |
| **Selected Features (Top-k)**   |198    | 1755    | 32 | 21 |

### evaluation of model sensitivity accross taxa

### visualise global and species specific FP/FN distribution
### analyse false positives
### analyse false negatives 



