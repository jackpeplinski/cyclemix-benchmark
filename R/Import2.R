require("SingleCellExperiment")
require("scater")
require("CycleMix")

counts <- read.table("./jackData/E-MTAB-2805.processed.1/G1_singlecells_counts.txt", header=TRUE)
rownames(counts) <- counts$EnsemblGeneID 
counts <- head(counts, -97) 
test <- DataFrame(row.names = counts$EnsemblGeneID, feature_symbol = counts$AssociatedGeneName)
sce <- SingleCellExperiment(assays=list(counts=counts),
                            rowData = DataFrame(row.names = counts$EnsemblGeneID,
                            feature_symbol = counts$AssociatedGeneName)) # is having not relevant colnames ok?

# Couldn't find the data set from the source 

# ignore
# ensemblGeneID <- counts$EnsemblGeneID
# genes <- counts$AssociatedGeneName
# counts <- subset(counts, select = -c(`EnsemblGeneID`,`EnsemblTranscriptID`, `AssociatedGeneName`, `GeneLength`) )