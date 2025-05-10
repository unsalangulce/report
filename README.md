# report
#   Single-Cell RNA-seq Analysis Report

**Dataset:** SCP2745_high_conf_CAS_cell_types.h5ad
**Date of Analysis:** (Generated from Notebook - Replace with actual date)
**Analysis Tool:** Scanpy, scvi-tools (via Python)

##  1.   Introduction

This report describes the initial steps in the analysis of a single-cell RNA sequencing (scRNA-seq) dataset, focusing on data loading, quality control, and doublet removal.

##  2.   Setup and Data Loading

###   2.1.  Library Import and Settings

The analysis was conducted using the Scanpy (1.10.4), and scvi-tools libraries in Python.

**Terminal Output:**

Libraries imported and Scanpy settings configured.

scanpy==1.10.4 anndata==0.11.3 umap==0.5.7 numpy==2.1.3 scipy==1.15.1 pandas==2.2.3 scikit-learn==1.6.1 statsmodels==0.14.4 python-igraph==NA leidenalg==NA


*(Note: I've filled in the Scanpy version from the notebook's output. You'll need to get the exact versions for the other libraries if needed using `print_header()` from Scanpy)*

###   2.2.  Data Loading

The dataset was loaded from the H5AD file: SCP2745_high_conf_CAS_cell_types.h5ad

**Terminal Output:**

Successfully loaded: SCP2745_high_conf_CAS_cell_types.h5ad

AnnData object summary:

AnnData object with n_obs × n_vars = 4000 × 33538
obs: 'cluster_set', 'high_conf_CAS_set_1_cell_types', 'high_conf_CAS_set_1_probabilities', 'high_conf_CAS_set_2_cell_types', 'high_conf_CAS_set_2_probabilities', 'high_conf_CAS_set_3_cell_types', 'high_conf_CAS_set_3_probabilities', 'granular_CAS_set_1_cell_types', 'granular_CAS_set_1_probabilities', 'granular_CAS_set_2_cell_types', 'granular_CAS_set_2_probabilities', 'granular_CAS_set_3_cell_types', 'granular_CAS_set_3_probabilities', 'biosample_id', 'donor_id', 'species', 'species__ontology_label', 'disease', 'disease__ontology_label', 'organ', 'organ__ontology_label', 'library_preparation_protocol', 'library_preparation_protocol__ontology_label', 'sex', 'total_mrna_umis'
var: 'gene_ids'
uns: 'cas_metadata'
obsm: 'X_pca', 'X_umap', 'cas_cl_scores'


The dataset comprises 4000 cells and 33538 genes/features.

###   2.3.  Initial Data Matrix Inspection

The data matrix `adata.X` has a shape of (4000, 33538).

**Terminal Output:**

(4000, 33538)


##  3.   Quality Control (QC) and Filtering

###   3.1.  Gene Filtering

Genes expressed in fewer than 10 cells were filtered out.

**Terminal Output:**

(Output of sc.pp.filter_genes not directly captured, but adata is modified in place)


*(Note: The exact number of removed genes isn't printed by `sc.pp.filter_genes` in the notebook, but the subsequent steps use the filtered `adata`)*

###   3.2.  Highly Variable Gene Selection

2000 highly variable genes were selected using the `seurat_v3` method.

**Terminal Output:**

(Output of sc.pp.highly_variable_genes not directly captured, but adata.var is modified in place)


*(Note: Again, the direct output of the function isn't explicitly printed, but the notebook proceeds using this selection.)*

##  4.   Doublet Removal

Doublets were predicted and removed using scvi-tools's SOLO method.

###   4.1.  scVI Model Training

An scVI model was trained on the filtered data.

**Terminal Output:**

Epoch 400/400: 100%|██████████| 400/400 [07:02<00:00, 1.06s/it, v_num=1, train_loss_step=827, train_loss_epoch=888]
Trainer.fit stopped: max_epochs=400 reached.


The scVI model was trained for 400 epochs.

###   4.2.  SOLO Model Training and Prediction

A SOLO model was trained to predict doublets, and then predictions were made.

**Terminal Output:**

Epoch 301/400: 75%|███████▌ | 301/400 [01:37<00:32, 3.08it/s, v_num=1, train_loss_step=0.198, train_loss_epoch=0.15]
Monitored metric validation_loss did not improve in the last 30 records. Best score: 0.152. Signaling Trainer to stop.


The SOLO model training stopped early due to lack of improvement in validation loss.

###   4.3.  Doublet Classification

The SOLO model classified cells as singlets or doublets.

**Terminal Output:**

                   doublet   singlet prediction
ACAGAAAAGGGAGGAC  5.365032e-02  0.946350    singlet
GGGTTATGTGGAACAC  6.602147e-02  0.933978    singlet
TGAGGAGTCCATTTCA  1.380594e-02  0.986194    singlet
CCTTTGGAGGAGTATT  9.305502e-01  0.069450    doublet
GAATCGTGTCAACATC  3.640850e-07  1.000000    singlet
...                        ...       ...        ...
GTAACACTCCTGGGTG  9.088891e-02  0.909111    singlet
GAGATGGTCCCAGCGA  5.978004e-02  0.940220    singlet
AGATCCAAGTTTAGGA  2.914435e-02  0.970856    singlet
TGAGTCAAGGTGGTTG  3.333474e-03  0.996666    singlet
CGTAAGTAGCTAGATA  9.774153e-01  0.022585    doublet

[4000 rows x 3 columns]


###   4.4.  Doublet Summary

A total of 327 cells were predicted to be doublets, while 3673 were classified as singlets.

**Terminal Output:**

        doublet  singlet
prediction
doublet         327      327
singlet        3673     3673


###   4.5.  Impact of Doublet Filtering

The doublet filtering step reduced the number of cells from 4000 to 3673, focusing subsequent analysis on high-quality singlets. This represents a removal of 327 cells identified as doublets.

**Terminal Output:**

        doublet  singlet
prediction
doublet         327      327
singlet        3673     3673


##  5.   Dimensionality Reduction

###   5.1.  Principal Component Analysis (PCA)

PCA was performed on the filtered data to reduce dimensionality.

**Terminal Output:**

(Output of sc.tl.pca not directly captured, but adata.obsm['X_pca'] is created)


*(Note: The notebook doesn't explicitly print PCA parameters, but `sc.tl.pca` is used, storing results in `adata.obsm['X_pca']`)*

###   5.2.  UMAP Embedding

UMAP embedding was computed for visualization of the cells in a 2-dimensional space.

**Terminal Output:**

computing UMAP
finished: added
'X_umap', UMAP coordinates (adata.obsm) (0:00:07)


UMAP coordinates were computed and stored in `adata.obsm['X_umap']`.

##  6.   Neighborhood Graph and Clustering

###   6.1.  Neighborhood Graph

A neighborhood graph was computed based on the PCA representation of the cells.

**Terminal Output:**

(Output of sc.pp.neighbors not directly captured, but adata.uns['neighbors'] is created)


###   6.2.  Leiden Clustering

Cells were clustered using the Leiden algorithm, based on the neighborhood graph.

**Terminal Output:**

running Leiden clustering
finished: found 13 clusters and added
'leiden', the cluster labels (adata.obs, categorical) (0:00:01)   


Leiden clustering identified 13 clusters.

##  7.   Visualization

###   7.1.  UMAP Visualization of Clusters

The clustered data was visualized on a UMAP plot, where each point represents a cell and is colored according to its Leiden cluster assignment.

*(Note: I cannot include the actual plot here. The notebook uses `sc.pl.umap(adata, color='leiden', legend_loc='on data')` to generate this plot. The report should reference the location of this figure within the notebook.)*

##  8.   Analysis: Initial Marker Gene Identification

The notebook begins the analysis phase by looking for marker genes that are characteristic of each cluster.

###   8.1.  Calculating doublet score difference

The notebook calculates the difference between doublet and singlet scores, and adds it to `adata.obs` as `dif`.

**Terminal Output:**

(Output of df['dif'] = df.doublet - df.singlet not directly printed, but adata.obs is modified)


###   8.2.  Visualization of doublet score difference

A plot is generated to show the distribution of the doublet score difference.

*(Note: The plot is generated using seaborn and matplotlib. Again, I can't include the plot itself, but the report should reference its location.)*

###   8.3.  Marker Gene Identification (Initial)

Initial steps towards marker gene identification were taken using `sc.tl.rank_genes_groups`.

**Terminal Output:**

(Output of sc.tl.rank_genes_groups not fully shown, but adata.uns['rank_genes_groups'] is created)


*(Note: The full output of `rank_genes_groups` is often extensive. The report should mention that this function was used and that the results are stored in `adata.uns['rank_genes_groups']`.)*
Note:
