# DFFB_persister_paper
Code for DFFB initial publication

This is the code for the initial DFFB persister cell publication. Code is for the single cell RNA sequencing analysis. 

Information about experiments:

A375 cells were cultured with 250 nM Dabrafenib and 25 nM Trametinib for 2 weeks to establish persister cells. PC9 cells were treated with 2.5 uM Erlotinib for 2 weeks. For both conditions, an untreated parental control was cultured alongside the drug-treated conditions. At the end of treatment, cells were lifted with trypsin and loaded onto a 10X Chromium instrument (10X Genomics) following the established protocol. Libraries were generated using the 10X Chromium Single Cell 3’ v3 kit as recommended. Quality control of the libraries was conducted with [whatever high sensitivity DNA kit] and then sequenced using NovaSeq S4. 

Read Alignment and data processing
Fastq files were aligned to the human “refdata-cellranger-GRCh38-3.0.0” genome with Cell Ranger version 3.1.0 with the “cellranger count” command to generate single cell feature counts for each library. The “filtered_feature_bc_matrix” generated for each population was used to create a “Seurat object” in the Seurat R package version 4.0.3. Cells containing greater than 1,000 and less then 7,500 features, and with less than 20% mitochondrial reads were included in downstream analyses. A cell cycle score was calculated for each cell using the default Seurat method, and this score was used to regress cell cycle during normalization and scaling with the SCTransform command. 

Determining variable features and mapping
The commands “RunPCA,” “RunUMAP,” “FindNeighbors,” and “FindClusters” were performed with default settings, with 30 dimensions used for “RunUMAP” and “FindNeighbors." In the “FindClusters” command, the resolution for A375 WT, A375 DFFB KO, and PC9 were set to 1.0, 0.4, and 0.2, respectfully, based on visualization of graphed clusters. The Seurat command “FindMarkerGenes” command was used with default parameters to identify differentially expressed genes between specified populations or clusters of cells. 

Gene-set Enrichment Analysis (GSEA)
GSEA was conducted with the ClusterProfiler R package (version 3.18.0) which calculated a normalized enrichment score (NES) using default parameters for each gene set. For the overlapping persister genes between A375 and PC9, genes differentially expressed in the same direction in both cell lines were analyzed with the GSEA/MSigDB website (http://www.gsea-msigdb.org/gsea) and a p-value and false discovery rate (FDR) was calculated with a hypergeometric test. The AUCell package version 1.12.0 was used to calculate gene set scores per cell, and the indicated thresholds were then visualized on a UMAP with Seurat. 

Individual gene expression analysis
The original count matrices were merged and filtered as previously described using Seurat. Expression values were then log-normalized and fold-change was calculated for a comparison between populations.
