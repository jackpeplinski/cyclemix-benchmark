library(Seurat)

exp.mat <- read.table(
    file = "./jackData/seurat/nestorawa_forcellcycle_expressionMatrix.txt",
    header = TRUE,
    as.is = TRUE, row.names = 1
)

print(head(exp.mat))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

marrow <- CreateSeuratObject(counts = exp.mat)
marrow <- NormalizeData(marrow)
marrow <- FindVariableFeatures(marrow, selection.method = "vst")
marrow <- ScaleData(marrow, features = rownames(marrow))
marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 6:10, nfeatures.print = 10)

marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
