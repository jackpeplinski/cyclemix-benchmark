convert_gene_symbols_to_ensembl_ids <- function(gene_symbols) {
    ensembl_ids_df <- suppressMessages(suppressWarnings(
        AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, columns = "ENSEMBL", keytype = "SYMBOL")
    ))

    ensembl_ids_df <- ensembl_ids_df[!is.na(ensembl_ids_df$ENSEMBL), ]

    # If there are multiple ENSEMBLs for a single SYMBOL, keep the first one
    ensembl_ids_df <- ensembl_ids_df[!duplicated(ensembl_ids_df$SYMBOL), ]

    return(ensembl_ids_df$ENSEMBL)
}

convert_ensemble_ids_to_gene_symbols <- function(ensembl_ids) {
    gene_symbols_df <- suppressMessages(suppressWarnings(
        AnnotationDbi::select(org.Hs.eg.db, keys = ensembl_ids, columns = "SYMBOL", keytype = "ENSEMBL")
    ))

    # If there are multiple SYMBOLs for a single ENSEMBL, keep the first one
    gene_symbols_df <- gene_symbols_df[!duplicated(gene_symbols_df$ENSEMBL), ]

    gene_symbols <- setNames(gene_symbols_df$SYMBOL, gene_symbols_df$ENSEMBL)

    gene_symbols <- gene_symbols[ensembl_ids]

    return(gene_symbols)
}
