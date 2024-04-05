library("biomaRt")
format_seurat_cc_to_mmus_ortho <- function() {
    data <- readRDS("./benchmarkData/SeuratCC_toMmus_ortho.rds")

    phase <- c("S")

    # Determine the number of rows in the data frame
    num_rows <- length(data$mmus_s)

    # Repeat the phase values to match the number of rows
    phase <- rep_len(phase, num_rows)

    # Combine the columns using cbind()
    result <- cbind(data$mmus_s, phase)

    # Convert the result to a data frame
    result <- as.data.frame(result)

    # Assign column names
    colnames(result) <- c("Gene", "Stage")

    ##### Phase 2
    phase <- c("G2M")

    # Determine the number of rows in the data frame
    num_rows <- length(data$mmus_g2m)

    # Repeat the phase values to match the number of rows
    phase <- rep_len(phase, num_rows)

    # Combine the columns using cbind()
    result2 <- cbind(data$mmus_g2m, phase)

    # Convert the result to a data frame
    result2 <- as.data.frame(result2)

    # Assign column names
    colnames(result2) <- c("Gene", "Stage")

    results <- rbind(result, result2)
    results$Gene <- results$Gene
    results$Dir <- 1
    results$Gene <- as.factor(results$Gene)
    results$Stage <- as.factor(results$Stage)

    ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

    ensembl_ids <- getBM(
        attributes = c("ensembl_gene_id", "external_gene_name"),
        filters = "external_gene_name",
        values = results$Gene,
        mart = ensembl
    )
    ensembl_ids$external_gene_name <- ensembl_ids$external_gene_name

    merged_df <- merge(results, ensembl_ids, by.x = "Gene", by.y = "external_gene_name", all.x = TRUE)

    rownames(merged_df) <- merged_df$ensembl_gene_id

    merged_df <- subset(merged_df, select = -ensembl_gene_id)
    return(merged_df)
}
MSeuratGeneSet <- format_seurat_cc_to_mmus_ortho()
