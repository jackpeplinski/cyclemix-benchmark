require("SingleCellExperiment")
require("scater")
require("CycleMix")
require("stringr")
library("Seurat")

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

    # remove unneeded info. counts is empty with only rownames = ensembl id now
    counts <- subset(counts,
        select = -c(
            id,
            gene.symbol,
            fpkm
        )
    )

    # duplicate counts variable for logcounts
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

    # for each file
    for (file in files) {

        # if the file name matches the col_data
        gsm <- substr(file, 1, 10)
        if (is.element(gsm, col_data_xml$gsm)) {
            # get the data for the file
            counts_data <- read.table(
                str_interp("./jackData/GSE42268/GSE42268_RAW/${file}"),
                header = TRUE
            )

            # format the data
            format_counts <- function(count_type) {

                # round data if counts
                if (count_type == "counts") {
                    data <- data.frame(
                        lapply(
                            counts_data,
                            function(x) if (is.numeric(x)) round(x) else x
                        )
                    )
                } else {
                    data <- counts_data
                }

                # sets rownames to ids
                rownames(data) <- data$id

                # remove uneeded data
                data <- subset(data,
                    select = -c(
                        id,
                        gene.symbol
                    )
                )

                colnames(data)[1] <- gsm
                return(data)
            }

            # the values for the gene expression
            logcounts_data <- format_counts("logcounts")
            counts_data <- format_counts("counts")
            # print(head(counts_data))
            # print(head(logcount_data))

            # the phases of cells
            col_data <- DataFrame(
                row.names = col_data_xml$gsm,
                Species = factor("Mus musculus"),
                cell_type1 = col_data_xml$cell_type1
            )
            # verify that the correct phase is matched
            # print(col_data)
            # print(col_data["GSM1036531", ])

            # combine phases and expressions
            logcounts <- cbind(logcounts, logcounts_data)
            counts <- cbind(counts, counts_data)
            # print(head(logcounts))
            # print(head(counts))
        }
    }

    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = col_data,
        rowData = row_data
    )
    # print(assays(sce))
    counts <- assay(sce, "counts")
    # libsizes <- colSums(logcounts)
    # size.factors <- libsizes / mean(libsizes)
    logcounts(sce) <- log2(t(t(logcounts)) + 1)
    # print(head(logcounts))
    # print(head(counts))
    # logcounts(sce) <- logcounts
    return(sce)
}

gse_output <- gse()
# head(logcounts(gse_output))
gse_classified <- classifyCells(gse_output, MGeneSets$Cyclone)
# gse_classified <- classifyCells(gse_output, subset(MGeneSets$Cyclone, Dir == 1)) # colData(gse_output)[50,] to view df with G2/M
summary(factor(gse_classified$phase))
table(factor(gse_classified$phase), gse_output$cell_type1)
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
emtab_2805_output <- emtab_2805() # should get

# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# emtab_2805_seurat <- as.Seurat(emtab_2805_output)
# emtab_2805_seurat <- NormalizeData(emtab_2805_seurat)
# emtab_2805_seurat <- FindVariableFeatures(emtab_2805_seurat, selection.method = "vst")
# emtab_2805_seurat <- ScaleData(emtab_2805_seurat, features = rownames(emtab_2805_seurat))
# emtab_2805_seurat <- RunPCA(emtab_2805_seurat, features = VariableFeatures(emtab_2805_seurat), ndims.print = 6:10, nfeatures.print = 10)
# CellCycleScoring(emtab_2805_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# look at changing the ensemble gene names to the gene symbols

# emtab_2805_classified <- classifyCells(emtab_2805_output, MGeneSets$Cyclone)
# summary(factor(emtab_2805_classified$phase))
# table(factor(emtab_2805_classified$phase), emtab_2805_output$cell_type1)
# plotMixture(emtab_2805_classified$fit[["G2M"]], BIC = TRUE)

cell_validity <- function() {
    invalid_row_counts <- c()
    for (i in colnames(logcounts(gse_output))) { # names of files
        # print(colnames(logcounts(gse_output)))
        invalid_row_count <- 0
        # if (i == "GSM1036501") {
        #     print(logcounts(gse_output)[i])
        # }
        file_log_counts <- logcounts(gse_output)[i]
        # print(colnames(fileLogCounts))
        colnames(file_log_counts)

        logcounts_rows <- row.names(file_log_counts)
        cyclone_rows <- row.names(MGeneSets$Cyclone)
        for (cycloneRow in cyclone_rows) {
            if (is.element(cycloneRow, logcounts_rows)) {
                # print(head(fileLogCounts)) # nolint
                if (file_log_counts[cycloneRow, ] < 5) {
                    invalid_row_count <- invalid_row_count + 1
                    # print(invalidRowCount) # nolint
                }
                # print(cycloneRow) # nolint
                # print(type(fileLogCounts[fileLogCounts[[1]] > 5, ])) # nolint
                # print(nrow(fileLogCounts)) # nolint
            }
        }
        invalid_row_counts <- append(invalid_row_counts, invalid_row_count
        / nrow(MGeneSets$Cyclone) * 100)
        # print(invalidRowCounts) # nolint
    }
    return(invalid_row_counts)
    # try subsetting MGeneSets$Cyclone by removing rows with Dir of -1 => this did improve it... but not a lot # nolint # nolint
    # look at average logcounts, the rows that have genes from MGeneSets$Cyclone. If the values are small (less than 5) they are not good indicators of the phase and how many are there. => 14 cells had poor logcounts, but 46 were not classified. # nolint
}
# validity <- cell_validity()

wilcox <- function() {
    "
    get row names from file
    for each ensemble id
        for each file
            get the value of that row name
            if that file is G2M store that value in a vector
            else store that value in another vector
    run a wilcox test
    "

    counts <- read.table(
        str_interp("./jackData/GSE42268/GSE42268_RAW/GSM1036480_EB5K_01.txt"),
        header = TRUE
    )
    ensemble_ids <- counts$id

    files <- list.files(
        path =
            "./jackData/GSE42268/GSE42268_RAW",
        pattern = ".*.txt"
    )

    col_data_xml <- data.frame(
        cell_type1 = c("G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "S", "S", "S", "S", "S", "S", "S", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1"), # nolint
        gsm = c("GSM1036483", "GSM1036484", "GSM1036485", "GSM1036486", "GSM1036487", "GSM1036488", "GSM1036489", "GSM1036490", "GSM1036491", "GSM1036492", "GSM1036493", "GSM1036494", "GSM1036495", "GSM1036496", "GSM1036500", "GSM1036501", "GSM1036502", "GSM1036503", "GSM1036504", "GSM1036505", "GSM1036506", "GSM1036507", "GSM1036508", "GSM1036509", "GSM1036510", "GSM1036511", "GSM1036512", "GSM1036513", "GSM1036522", "GSM1036523", "GSM1036524", "GSM1036525", "GSM1036526", "GSM1036527", "GSM1036528", "GSM1036529", "GSM1036530", "GSM1036531", "GSM1036532", "GSM1036533", "GSM1036534", "GSM1036535", "GSM1036536", "GSM1036537", "GSM1036538", "GSM1036539", "GSM1036540", "GSM1036541", "GSM1036542", "GSM1036543", "GSM1036544", "GSM1036545", "GSM1036546", "GSM1036547", "GSM1036548", "GSM1036549", "GSM1036550", "GSM1036551", "GSM1036552", "GSM1036553", "GSM1036554", "GSM1036555", "GSM1036556") # nolint
    )

    for (ensemble_id in ensemble_ids) {
        g2m <- c()
        g1 <- c()
        s <- c()
        for (file in files) {
            gsm <- substr(file, 1, 10)
            counts <- read.table(
                str_interp("./jackData/GSE42268/GSE42268_RAW/${file}"),
                header = TRUE
            )
            fpkm <- counts[counts$id == ensemble_id, ]$fpkm
            # print(col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1)
            # print(length(col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1) > 0)
            if (length(col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1) > 0 &&
                col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1 == "G2/M") {
                g2m <- append(g2m, fpkm)
            }
            if (length(col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1) > 0 &&
                col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1 == "S") {
                g2m <- append(s, fpkm)
            } else {
                g1 <- append(g1, fpkm)
            }
        }
        print(paste(ensemble_id, wilcox.test(g2m, s)$p.value))
        print(paste(ensemble_id, wilcox.test(g2m, g1)$p.value))
    }
}
# wilcox()
# for each gene in the cyclemix
# look in the apply function
# wilcox.test(#vector of all the values of the cells that are G2M, #same but for S/G1) # nolint

wilcox2 <- function() {
    "
            row.name (ensembleID)    gsmcellname
            ENSMUSG00000000049       [fpkm]
            phase                    [phase]
    "

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
    ensemble_ids <- counts$id

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
    # logcounts <- counts

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

        # if file name and col data are the same
        if (is.element(gsm, col_data_xml$gsm)) {
            counts_data <- read.table(
                str_interp("./jackData/GSE42268/GSE42268_RAW/${file}"),
                header = TRUE
            )

            format_counts <- function(count_type) {
                # if (count_type != "logcounts") {
                #     data <- data.frame(
                #         lapply(
                #             counts_data,
                #             function(x) if (is.numeric(x)) round(x) else x
                #         )
                #     )
                # } else {
                data <- counts_data
                # }

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
            #         logcounts_data <- format_counts("logcounts")
            counts_data <- format_counts("counts")

            col_data <- DataFrame(
                row.names = col_data_xml$gsm,
                Species = factor("Mus musculus"),
                cell_type1 = col_data_xml$cell_type1 # why the "1" at the end?
            )
            # verify that the correct phase is matched
            # print(col_data)
            # print(col_data["GSM1036531", ])

            #         logcounts <- cbind(logcounts, logcounts_data)
            counts <- cbind(counts, counts_data)
        }
    }
    counts[nrow(counts) + 1, ] <- list("G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "S", "S", "S", "S", "S", "S", "S", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1")
    rownames(counts)[36808] <- "phase"
    names(df)[which(df == 1, arr.ind = T)[, "col"]]

    # for (ensemble_id in ensemble_ids) {
    #     g2m <- c()
    #     g1 <- c()
    #     s <- c()

    #     g2m <- counts

    #     print(paste(ensemble_id, wilcox.test(g2m, s)$p.value))
    #     print(paste(ensemble_id, wilcox.test(g2m, g1)$p.value))
    # }
    # sce <- SingleCellExperiment(
    #     assays = list(counts = counts),
    #     colData = col_data, # look at split function
    #     rowData = row_data
    # )
    # # print(assays(sce))
    # counts <- assay(sce, "counts") # double check this and the logcounts is what she wanted. # nolint
    # # print(head(logcounts))
    # logcounts(sce) <- logcounts
    # return(sce)
}
# wilcox2()
