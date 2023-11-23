source("benchmark.R")

test_synthesize_gse_42268 <- function(sce) {
    c1 <- c(1, 2, 3)
    c2 <- c(4, 5, 6)
    c3 <- c(7, 8, 9)
    c <- c(c1, c2, c3)
    counts <- matrix(c, ncol = 3, nrow = 3)
    gene_names <- c("a", "b", "c")
    rownames(counts) <- gene_names
    cell_names <- c("cell1", "cell2", "cell3")
    colnames(counts) <- cell_names

    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        rowData = DataFrame(row.names = gene_names, feature_symbol = factor(gene_names)),
        colData = DataFrame(row.names = cell_names, cell_type1 = c(factor("G1"), factor("G1"), factor("G1")))
    )
    actual <- synthesize_gse_42268_by_type(sce, "G1")

    expected <- matrix(c(2.5, 3.5, 4.5, 4, 5, 6, 5.5, 6.5, 7.5, c), ncol = 6, nrow = 3)
    rownames(expected) <- gene_names
    col_names <- c("cell1-2", "cell1-3", "cell2-3", cell_names)
    colnames(expected) <- col_names
    return(actual == expected)
}
