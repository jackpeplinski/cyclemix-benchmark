require("SingleCellExperiment")
require("scater")
require("CycleMix")
require("stringr")

# Don't need to log normalize, each file is a cell.
# Round fpkm for counts, and put fpkm into the logcounts
gse <- function() {
    # counts should be something like
    # with headers: row.name (ensembleID),   gsmcellname
    #               e.g.,ENSMUSG00000000049  [fpkm]

    # load counts info; doesn't matter which file, could be any
    counts <- read.table(
        str_interp("./jackData/GSE42268/GSE42268_RAW/GSM1036480_EB5K_01.txt"),
        header = TRUE
    )

    # adds row.name with ensembleID
    rownames(counts) <- counts$id

    # adds gene.symbol with gene symbol
    row_data <- DataFrame(
        row.names = counts$id,
        feature_symbol = factor(counts$gene.symbol)
    )

    # remove unneeded info
    counts <- subset(counts,
        select = -c(
            id,
            gene.symbol,
            fpkm
        )
    )
    logcounts <- counts

    # copy/pasted this from gseParse
    col_data_xml <- data.frame(
        cell_type1 = c("G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "S", "S", "S", "S", "S", "S", "S", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1"), # nolint
        gsm = c("GSM1036483", "GSM1036484", "GSM1036485", "GSM1036486", "GSM1036487", "GSM1036488", "GSM1036489", "GSM1036490", "GSM1036491", "GSM1036492", "GSM1036493", "GSM1036494", "GSM1036495", "GSM1036496", "GSM1036500", "GSM1036501", "GSM1036502", "GSM1036503", "GSM1036504", "GSM1036505", "GSM1036506", "GSM1036507", "GSM1036508", "GSM1036509", "GSM1036510", "GSM1036511", "GSM1036512", "GSM1036513", "GSM1036522", "GSM1036523", "GSM1036524", "GSM1036525", "GSM1036526", "GSM1036527", "GSM1036528", "GSM1036529", "GSM1036530", "GSM1036531", "GSM1036532", "GSM1036533", "GSM1036534", "GSM1036535", "GSM1036536", "GSM1036537", "GSM1036538", "GSM1036539", "GSM1036540", "GSM1036541", "GSM1036542", "GSM1036543", "GSM1036544", "GSM1036545", "GSM1036546", "GSM1036547", "GSM1036548", "GSM1036549", "GSM1036550", "GSM1036551", "GSM1036552", "GSM1036553", "GSM1036554", "GSM1036555", "GSM1036556") # nolint
    )

    files <- list.files(
        path =
            "./jackData/GSE42268/GSE42268_RAW",
        pattern = ".*.txt"
    )

    for (file in files) {
        gsm <- substr(file, 1, 10)
        if (is.element(gsm, col_data_xml$gsm)) {
            counts_data <- read.table(
                str_interp("./jackData/GSE42268/GSE42268_RAW/${file}"),
                header = TRUE
            )

            format_counts <- function(count_type) {
                if (count_type != "logcounts") {
                    data <- data.frame(
                        lapply(
                            counts_data,
                            function(x) if (is.numeric(x)) round(x) else x
                        )
                    )
                } else {
                    data <- counts_data
                }

                rownames(data) <- data$id

                data <- subset(data,
                    select = -c(
                        id,
                        gene.symbol
                    )
                )

                colnames(data)[1] <- gsm
                return(data)
            }
            logcounts_data <- format_counts("logcounts")
            counts_data <- format_counts("counts")

            col_data <- DataFrame(
                row.names = col_data_xml$gsm,
                Species = factor("Mus musculus"),
                cell_type1 = col_data_xml$cell_type1 # why the "1" at the end?
                # Source = factor("ESC") # this may be different
            )
            # verify that the correct phase is matched
            # print(col_data)
            # print(col_data["GSM1036531", ])

            logcounts <- cbind(logcounts, logcounts_data)
            counts <- cbind(counts, counts_data)
        }
    }

    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = col_data, # look at split function
        rowData = row_data
    )
    print(assays(sce))
    counts <- assay(sce, "counts") # double check this and the logcounts is what she wanted.
    # print(head(logcounts))
    logcounts(sce) <- logcounts
    return(sce)
}

gse_output <- gse() # head(logcounts(gse_output))
# gse_classified <- classifyCells(gse_output, MGeneSets$Cyclone)
# gse_classified <- classifyCells(gse_output, subset(MGeneSets$Cyclone, Dir == 1)) # colData(gse_output)[50,] to view df with G2/M
# summary(factor(gse_classified$phase))
# table(factor(gse_classified$phase), gse_output$cell_type1)
# plotMixture(gse_classified$fit[["G2M"]], BIC = TRUE)
# plotMixture(gse_classified$fit[["S"]], BIC = TRUE)
# plotMixture(gse_classified$fit[["G1"]], BIC = TRUE)

emtab_2805 <- function() {
    emtab_2805_file <- function(file_name) {
        counts <- read.table(
            str_interp("./jackData/E-MTAB-2805/E-MTAB-2805.processed.1/${file_name}.txt"), # nolint
            header = TRUE
        )

        rownames(counts) <- counts$EnsemblGeneID

        counts <- head(counts, -97)
        row_data <- DataFrame(
            row.names = counts$EnsemblGeneID,
            feature_symbol = factor(counts$AssociatedGeneName)
        )

        col_data <- DataFrame(
            row.names = colnames(counts)[5:100],
            Species = factor("Mus musculus"),
            cell_type1 = factor(
                substr(colnames(counts)[5:100], 1, 2)
            ),
            Source = factor("ESC")
        )

        counts <- subset(counts,
            select = -c(
                EnsemblGeneID,
                EnsemblTranscriptID,
                AssociatedGeneName,
                GeneLength
            )
        )
        # print(head(counts))
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
# emtab_2805_output <- emtab_2805() # should get
# emtab_2805_classified <- classifyCells(emtab_2805_output, MGeneSets$Cyclone)
# summary(factor(emtab_2805_classified$phase))
# table(factor(emtab_2805_classified$phase), emtab_2805_output$cell_type1)
# plotMixture(emtab_2805_classified$fit[["G2M"]], BIC = TRUE)

cellValidity <- function() {
    invalidRowCounts <- c()
    for (i in colnames(logcounts(gse_output))) {
        invalidRowCount <- 0
        # if (i == "GSM1036501") {
        #     print(logcounts(gse_output)[i])
        # }
        fileLogCounts <- logcounts(gse_output)[i]
        # print(colnames(fileLogCounts))
        colnames(fileLogCounts)

        logcountsRows <- row.names(fileLogCounts)
        cycloneRows <- row.names(MGeneSets$Cyclone)
        for (cycloneRow in cycloneRows) {
            if (is.element(cycloneRow, logcountsRows)) {
                # print(head(fileLogCounts))
                if (fileLogCounts[cycloneRow, ] < 5) {
                    invalidRowCount <- invalidRowCount + 1
                    # print(invalidRowCount)
                }
                # print(cycloneRow)
                # print(type(fileLogCounts[fileLogCounts[[1]] > 5, ]))
                # print(nrow(fileLogCounts))
            }
        }
        invalidRowCounts <- append(invalidRowCounts, invalidRowCount
        / nrow(MGeneSets$Cyclone) * 100)
        # print(invalidRowCounts)
    }
    return(invalidRowCounts)
    # try subsetting MGeneSets$Cyclone by removing rows with Dir of -1 => this did improve it... but not a lot
    # look at average logcounts, the rows that have genes from MGeneSets$Cyclone. If the values are small (less than 5) they are not good indicators of the phase and how many are there. => 14 cells had poor logcounts, but 46 were not classified.
}
validity <- cellValidity()

wilcox <- function() {
    # counts should be something like
    # with headers: row.name (ensembleID),   gsmcellname
    #               e.g.,ENSMUSG00000000049  [fpkm]

    # load counts info; doesn't matter which file, could be any
    counts <- read.table(
        str_interp("./jackData/GSE42268/GSE42268_RAW/GSM1036480_EB5K_01.txt"),
        header = TRUE
    )

    # adds row.name with ensembleID
    rownames(counts) <- counts$id

    # adds gene.symbol with gene symbol
    row_data <- DataFrame(
        row.names = counts$id,
        feature_symbol = factor(counts$gene.symbol)
    )

    # remove unneeded info
    counts <- subset(counts,
        select = -c(
            id,
            gene.symbol,
            fpkm
        )
    )
    logcounts <- counts

    # copy/pasted this from gseParse
    col_data_xml <- data.frame(
        cell_type1 = c("G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "S", "S", "S", "S", "S", "S", "S", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1"), # nolint
        gsm = c("GSM1036483", "GSM1036484", "GSM1036485", "GSM1036486", "GSM1036487", "GSM1036488", "GSM1036489", "GSM1036490", "GSM1036491", "GSM1036492", "GSM1036493", "GSM1036494", "GSM1036495", "GSM1036496", "GSM1036500", "GSM1036501", "GSM1036502", "GSM1036503", "GSM1036504", "GSM1036505", "GSM1036506", "GSM1036507", "GSM1036508", "GSM1036509", "GSM1036510", "GSM1036511", "GSM1036512", "GSM1036513", "GSM1036522", "GSM1036523", "GSM1036524", "GSM1036525", "GSM1036526", "GSM1036527", "GSM1036528", "GSM1036529", "GSM1036530", "GSM1036531", "GSM1036532", "GSM1036533", "GSM1036534", "GSM1036535", "GSM1036536", "GSM1036537", "GSM1036538", "GSM1036539", "GSM1036540", "GSM1036541", "GSM1036542", "GSM1036543", "GSM1036544", "GSM1036545", "GSM1036546", "GSM1036547", "GSM1036548", "GSM1036549", "GSM1036550", "GSM1036551", "GSM1036552", "GSM1036553", "GSM1036554", "GSM1036555", "GSM1036556") # nolint
    )

    files <- list.files(
        path =
            "./jackData/GSE42268/GSE42268_RAW",
        pattern = ".*.txt"
    )

    for (file in files) {
        gsm <- substr(file, 1, 10)
        if (is.element(gsm, col_data_xml$gsm)) {
            counts_data <- read.table(
                str_interp("./jackData/GSE42268/GSE42268_RAW/${file}"),
                header = TRUE
            )

            format_counts <- function(count_type) {
                if (count_type != "logcounts") {
                    data <- data.frame(
                        lapply(
                            counts_data,
                            function(x) if (is.numeric(x)) round(x) else x
                        )
                    )
                } else {
                    data <- counts_data
                }

                rownames(data) <- data$id

                data <- subset(data,
                    select = -c(
                        id,
                        gene.symbol
                    )
                )

                colnames(data)[1] <- gsm
                return(data)
            }
            logcounts_data <- format_counts("logcounts")
            counts_data <- format_counts("counts")

            col_data <- DataFrame(
                row.names = col_data_xml$gsm,
                Species = factor("Mus musculus"),
                cell_type1 = col_data_xml$cell_type1 # why the "1" at the end?
                # Source = factor("ESC") # this may be different
            )
            # verify that the correct phase is matched
            # print(col_data)
            # print(col_data["GSM1036531", ])

            logcounts <- cbind(logcounts, logcounts_data)
            counts <- cbind(counts, counts_data)
        }
    }

    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = col_data, # look at split function
        rowData = row_data
    )
    counts <- assay(sce, "counts") # double check this and the logcounts is what she wanted.
    # print(head(logcounts))
    logcounts(sce) <- logcounts
    return(sce)
}

# for each gene
# wilcox.test(#vector of all the values of the cells that are G2M , #same but for S/G1)
