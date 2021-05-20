#Load libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)

## Load the A375 parental dataset
wt.par.data <- Read10X(data.dir = "~/Desktop/apar_output_files/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
apar <- CreateSeuratObject(counts = wt.par.data, project = "1_A375_parental", min.cells = 3, min.features = 200)
# Load the A375 persister dataset
wt.pers.data <- Read10X(data.dir = "~/Desktop/apers_output_files/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
apers <- CreateSeuratObject(counts = wt.pers.data, project = "2_A375_persisters", min.cells = 3, min.features = 200)
# Initialize the Seurat object with the raw (non-normalized data).
all.combined <- merge(apar, y = c(apers), add.cell.ids = c("wt_par", "wt_pers"), project = "A375_par_pers_merged")
all.combined[["percent.mt"]] <- PercentageFeatureSet(object = all.combined, pattern = "^MT-")
all.combined <- subset(x = all.combined, subset = nFeature_RNA > 1000 & percent.mt < 20 & nFeature_RNA < 7500)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.combined <- CellCycleScoring(all.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
all.combined <- SCTransform(all.combined, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
all.combined <- RunPCA(all.combined, verbose = FALSE)
all.combined <- RunUMAP(all.combined, dims = 1:30, verbose = FALSE)
all.combined <- FindNeighbors(all.combined, dims = 1:30, verbose = FALSE)
all.combined <- FindClusters(all.combined, verbose = FALSE, resolution = 1.0)

#Create a UMAP
DimPlot(all.combined, label = FALSE, group.by = "orig.ident", label.size = 8)

#Extract expression matrix
a375avgexpr <- AverageExpression(all.combined, group.by = "orig.ident", return.seurat = TRUE)
exprMat <- GetAssayData(object = a375avgexpr, slot = "data")
write.csv(exprMat, "~/Desktop/2021.4.12_A375_par_per_dtep_matrix.csv")

#Get breakdown of cell cycle
parental.cells <- subset(all.combined, subset = (orig.ident == "1_A375_parental"))
table(parental.cells$Phase)

#Identify marker genes
par.markers <- FindMarkers(object = combined.sct, ident.1 = "1_A375_parental", ident.2 = "2_A375_persisters", group.by = "orig.ident", min.pct = 0.20)

#Code for ClusterProfiler
#Set up list of genes
d <- read.csv("/Users/augustwilliams/Desktop/2021.1.23_ko_cluster_0_v4_cells.csv")
geneList <- d[,3]
names(geneList) <- as.character(d[,1])
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)
gene <- names(geneList)[abs(geneList) > 0.301]
head(gene)

#Load GSEA gene sets
library(msigdbr)
msigdbr_species()
#Obtain a specific collection from msigdbr
#H: hallmark gene sets
#C1: positional gene sets
#C2: curated gene sets
#C3: motif gene sets
#C4: computational gene sets
#C5: GO gene sets
#C6: oncogenic signatures
#C7: immunologic signatures
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, human_gene_symbol)
go <- msigdbr(species = "Homo sapiens", category = "C5") %>% dplyr::select(gs_name, human_gene_symbol)
oncosigs <- msigdbr(species = "Homo sapiens", category = "C6") %>% dplyr::select(gs_name, human_gene_symbol)
head(m_t2g)

#Enter list and gene set into cluster profiler
library(clusterProfiler)
em2 <- GSEA(geneList, TERM2GENE = oncosigs, minGSSize =5)
em2
library(enrichplot)
library(ggplot2)
em2@result
data <- data.frame(
  name=c(em2$ID) ,  
  value=c(em2$NES)
)
#generate plot for GSEA normalized enrichment score
library(forcats)
data %>%
  mutate(name = fct_reorder(name, value)) %>%
  ggplot( aes(x=name, y=value)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()






