#AUCell pipeline
library(AUCell)
counts_matrix <- GetAssayData(object = all.combined, slot = "counts")
persister_rankings <- AUCell_buildRankings(counts_matrix, nCores=1, plotStats=TRUE)
save(persister_rankings, file="A375_persisters_cell_rankings.RData")
library(GSEABase)
gmtFile <- paste("/Users/augustwilliams/Desktop/h.all.v7.2.symbols.gmt")
geneSets <- getGmt(gmtFile)
geneSets <- subsetGeneSets(geneSets, rownames(counts_matrix))
cbind(nGenes(geneSets))
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep=""))
#The actual command
cells_AUC <- AUCell_calcAUC(geneSets, persister_rankings)
save(cells_AUC, file="A375_immunopers_cells_AUC.RData")
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)

#create assignment to plot on UMAP
apop <- cells_assignment$`HALLMARK_APOPTOSIS (156g)`$assignment

#Generate UMAP colored with assignment
DimPlot(object = combined.sct, reduction = 'umap', label = FALSE, cells.highlight = apop)
