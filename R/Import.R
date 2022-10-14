counts00 <- read.table("./jackData/E-MTAB-2805.processed.1/G1_singlecells_counts.txt", header=TRUE)
counts01 <- data.frame(row.names = counts00$EnsemblGeneID, Gene = counts00$AssociatedGeneName, Stage = "G1")

counts02 <- read.table("./jackData/E-MTAB-2805.processed.1/G2M_singlecells_counts.txt", header=TRUE)
counts03 <- data.frame(row.names = counts02$EnsemblGeneID, Gene = counts02$AssociatedGeneName, Stage = "G2M")

counts04 <- read.table("./jackData/E-MTAB-2805.processed.1/S_singlecells_counts.txt", header=TRUE)
counts05 <- data.frame(row.names = counts04$EnsemblGeneID, Gene = counts04$AssociatedGeneName, Stage = "S1")

test <- rbind(counts01, counts03, counts05)


