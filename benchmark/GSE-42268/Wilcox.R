source("./benchmark/GSE-42268/Format.R")
# for each gene in the GSE-42268 dataset
# wilcox.test(#vector of all the values of the cells that are G2M, #same but for S/G1)

get_wilcox_fp <- function(phase1, phase2) {
    return(paste0("./output/wilcox-", gsub("/", "", phase1), "-", gsub("/", "", phase2), ".txt"))
}

get_counts_by_phase <- function() {
    sce <- get_sce()
    phases <- c("S", "G1", "G2/M")
    counts_by_phase <- list()
    for (phase in phases) {
        cell_names_by_phase <- rownames(colData(sce)[colData(sce)$cell_type1 == phase, , drop = FALSE])
        counts_by_phase[[phase]] <- logcounts(sce)[, cell_names_by_phase]
    }
    return(counts_by_phase)
}

write_wilcox_output <- function(phase1, phase2, counts_by_phase, wilcox_fp) {
    counts_by_phase1 <- counts_by_phase[[phase1]]
    counts_by_phase2 <- counts_by_phase[[phase2]]
    print("Writing to file...")
    file <- file(wilcox_fp, "w")
    writeLines("id p_value", file)
    for (row in rownames(counts_by_phase1)) {
        x <- counts_by_phase1[row, ]
        y <- counts_by_phase2[row, ]
        # also display the mean and output to a file so that we can sort it by pvalue
        pval <- p.adjust(wilcox.test(as.numeric(unlist((x))), as.numeric(unlist(y)))$p.value, method = "fdr")
        line <- paste(row, pval, sep = " ")
        writeLines(line, file)
    }
    close(file)
    print("Writing done.")
}

read_wilcox_output <- function(wilcox_fp) {
    df <- read.table(
        wilcox_fp,
        header = TRUE
    )
    # how many are below 0.05
    print(paste0(nrow(df[df$p_value < 0.05, ]), " genes have a p-value below 0.05 out of ", nrow(df)))
    # print(paste0(mean(df$p_value), " is the mean p-value"))
}

counts_by_phase <- get_counts_by_phase()
phase1 <- "G2/M"
phase2 <- "S"
wilcox_fp <- get_wilcox_fp(phase1, phase2)
write_wilcox_output(phase1, phase2, counts_by_phase, wilcox_fp)
read_wilcox_output(wilcox_fp)
