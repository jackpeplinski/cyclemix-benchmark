require("SingleCellExperiment")
require("scater")
require("CycleMix")
require("stringr")
library("Seurat")
library("biomaRt")
library("org.Mm.eg.db")
options(max.print = 20)

format_emtab_2805 <- function() {
    emtab_2805_file <- function(file_name) {
        # get values
        counts <- read.table(
            str_interp("./benchmarkData/E-MTAB-2805/E-MTAB-2805.processed.1/${file_name}.txt"),
            header = TRUE
        )
        counts$AssociatedGeneName <- counts$AssociatedGeneName
        counts <- counts[!duplicated(counts$AssociatedGeneName), ]
        counts <- na.omit(counts)

        # set rownames
        rownames(counts) <- counts$AssociatedGeneName

        # remove unused values
        counts <- head(counts, -97)

        # build needed df for sce
        row_data <- DataFrame(
            row.names = counts$AssociatedGeneName,
            feature_symbol = factor(counts$AssociatedGeneName)
        )

        # build needed df for sce
        col_data <- DataFrame(
            row.names = colnames(counts)[5:100],
            Species = factor("Mus musculus"),
            cell_type1 = factor(
                substr(colnames(counts)[5:100], 1, 2)
            ),
            Source = factor("ESC")
        )

        # remove unneeded data
        counts <- subset(counts,
            select = -c(
                EnsemblGeneID,
                EnsemblTranscriptID,
                AssociatedGeneName,
                GeneLength
            )
        )

        # build sce
        sce <- SingleCellExperiment(
            assays = list(counts = counts),
            colData = col_data,
            rowData = row_data
        )
        counts <- assay(sce, "counts")
        libsizes <- colSums(counts)
        size.factors <- libsizes / mean(libsizes)
        logcounts(sce) <- log2(t(t(counts) / size.factors) + 1)
        return(sce)
    }

    emtab_2805_res <- cbind(
        emtab_2805_file("G1_singlecells_counts"),
        emtab_2805_file("G2M_singlecells_counts"),
        emtab_2805_file("S_singlecells_counts")
    )
    return(emtab_2805_res)
}

classify_emtab_2805 <- function() {
    cat("===EMTAB 2805 | CycleMix | MGeneSets$Cyclone===\n")
    emtab_sce <<- format_emtab_2805()
    emtab_cm_cy <<- classifyCells(emtab_sce, MGeneSets$Cyclone)
    print(table(factor(emtab_cm_cy$phase), emtab_sce$cell_type1))
    cat("===EMTAB 2805 | CycleMix | MSeuratGeneSet===\n")
    emtab_cm_se <<- classifyCells(emtab_sce, MSeuratGeneSet)
    # plotMixture(emtab_cm_se$fit[["G1"]], BIC = TRUE)
    # plotMixture(emtab_cm_se$fit[["G2M"]], BIC = TRUE)
    # plotMixture(emtab_cm_se$fit[["S"]], BIC = TRUE)
    emtab_cm_se_table <- table(factor(emtab_cm_se$phase), emtab_sce$cell_type1)
    print(emtab_cm_se_table)

    emtab_seurat <- as.Seurat(emtab_sce)
    emtab_seurat <- NormalizeData(emtab_seurat)
    emtab_seurat <- FindVariableFeatures(emtab_seurat, selection.method = "vst")
    emtab_seurat <- ScaleData(emtab_seurat, features = rownames(emtab_seurat))
    emtab_seurat <<- RunPCA(emtab_seurat, features = VariableFeatures(emtab_seurat), ndims.print = 6:10, nfeatures.print = 10)
    cat("===EMTAB 2805 | Seurat | MGeneSets$Cyclone===\n")
    s.genes <- MGeneSets$Cyclone$Gene[MGeneSets$Cyclone$Stage == "S"]
    g2m.genes <- MGeneSets$Cyclone$Gene[MGeneSets$Cyclone$Stage == "G2M"]
    emtab_seurat_cy <<- CellCycleScoring(emtab_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    print(table(emtab_seurat_cy[[]]$Phase, emtab_seurat[[]]$orig.ident))
    # plotMixture(emtab_cm_cy$fit[["G1"]], BIC = TRUE)
    # plotMixture(emtab_cm_c$fit[["G2M"]], BIC = TRUE)
    # plotMixture(emtab_cm_cy$fit[["S"]], BIC = TRUE)
    cat("===EMTAB 2805 | Seurat | MSeuratGeneSet===\n")
    s.genes <- seurat_mouse_orth$mmus_s
    g2m.genes <- seurat_mouse_orth$mmus_g2m
    emtab_seurat_se <<- CellCycleScoring(emtab_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    print(table(emtab_seurat_se[[]]$Phase, emtab_seurat[[]]$orig.ident))
}
classify_emtab_2805()
