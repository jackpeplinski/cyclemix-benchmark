source("./benchmark/GSE-42268/Format.r")
require("CycleMix")
library("Seurat")

MSeuratGeneSet <- readRDS("./benchmarkData/MSeuratGeneSet.RDS")
seurat_mouse_orth <- readRDS("./benchmarkData/SeuratCC_toMmus_ortho.rds")

classify_gse_42268 <- function(gse_sce) {
    cat("===GSE 42268 | CycleMix | MGeneSets$Cyclone===\n")
    gse_cm_cy <- classifyCells(gse_sce, MGeneSets$Cyclone)
    print(table(factor(gse_cm_cy$phase), gse_sce$cell_type1))

    cat("===GSE 42268 | CycleMix | MSeuratGeneSet===\n")
    gse_cm_se <- classifyCells(gse_sce, MSeuratGeneSet)
    print(table(factor(gse_cm_se$phase), gse_sce$cell_type1))

    cat("===Seurat===\n")
    gse_seurat <- as.Seurat(gse_sce)
    gse_seurat <- NormalizeData(gse_seurat)
    gse_seurat <- FindVariableFeatures(gse_seurat, selection.method = "vst")
    gse_seurat <- ScaleData(gse_seurat, features = rownames(gse_seurat))
    gse_seurat <- RunPCA(gse_seurat, features = VariableFeatures(gse_seurat), ndims.print = 6:10, nfeatures.print = 10)

    cat("===GSE 42268 | Seurat | MGeneSets$Cyclone===\n")
    s.genes <- MGeneSets$Cyclone$Gene[MGeneSets$Cyclone$Stage == "S"]
    s.genes <- convert_mgene_symbols_to_ensembl_ids(as.character(s.genes))
    g2m.genes <- MGeneSets$Cyclone$Gene[MGeneSets$Cyclone$Stage == "G2M"]
    g2m.genes <- convert_mgene_symbols_to_ensembl_ids(as.character(g2m.genes))
    gse_seurat_cy <- CellCycleScoring(gse_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    print(table(gse_seurat_cy[[]]$Phase, gse_seurat[[]]$cell_type1))

    cat("===GSE 42268 | Seurat | MSeuratGeneSet===\n")
    s.genes <- convert_mgene_symbols_to_ensembl_ids(seurat_mouse_orth$mmus_s)
    g2m.genes <- convert_mgene_symbols_to_ensembl_ids(seurat_mouse_orth$mmus_g2m)
    gse_seurat_se <- CellCycleScoring(gse_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    print(table(gse_seurat_se[[]]$Phase, gse_seurat[[]]$cell_type1))
}

gse_sce <- get_sce()
cat("***Non-synthetic***\n")
classify_gse_42268(gse_sce)
# cat("***Synthetic***\n")
# classify_gse_42268(synthesize_gse_42268(gse_sce))
