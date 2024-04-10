# for each gene in the cyclemix
# wilcox.test(#vector of all the values of the cells that are G2M, #same but for S/G1)
wilcox_gse_42268 <- function() {
    # this may be broken
    "
    get row names from file
    for each ensemble id
        for each file
            get the value of that row name
            if that file is G2M store that value in a vector
            else store that value in another vector
    run a wilcox test
    "

    # get ensemble_ids
    counts <- read.table(
        str_interp("./benchmarkData/GSE42268/GSE42268_RAW/GSM1036480_EB5K_01.txt"),
        header = TRUE
    )
    ensemble_ids <- counts$id

    # get list of files
    files <- list.files(
        path =
            "./benchmarkData/GSE42268/GSE42268_RAW",
        pattern = ".*.txt"
    )

    col_data_xml <- readRDS("./benchmarkData/ColDataXML.rds")


    # for each id
    for (ensemble_id in ensemble_ids) {
        # if the element is in Cyclone set
        if (is.element(ensemble_id, rownames(MGeneSets$Cyclone))) {
            g2m <- c()
            g1 <- c()
            s <- c()

            # for each file
            for (file in files) {
                # get the file name
                gsm <- substr(file, 1, 10)

                # get values for the file
                counts <- read.table(
                    str_interp("./benchmarkData/GSE42268/GSE42268_RAW/${file}"),
                    header = TRUE
                )

                # get the values for that ID
                fpkm <- counts[counts$id == ensemble_id, ]$fpkm

                # print(col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1)
                # print(length(col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1) > 0)

                # add the fpkm value for the ensemble for that cell
                if (length(col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1) > 0 &&
                    col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1 == "G2/M") {
                    g2m <- append(g2m, fpkm)
                } else if (length(col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1) > 0 &&
                    col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1 == "S") {
                    g2m <- append(s, fpkm)
                } else {
                    g1 <- append(g1, fpkm)
                }
            }
            print(g2m)
            print(paste(ensemble_id, wilcox.test(g2m, s)$p.value))
            print(paste(ensemble_id, wilcox.test(g2m, g1)$p.value))
        }
    }
}

wilcox_fast_gse_42268 <- function() {
    "
        row.name (ensembleID)    gsmcellname
        ENSMUSG00000000049       0.002000
    "

    # set row names to ensemble ids
    counts <- read.table(
        str_interp("./benchmarkData/GSE42268/GSE42268_RAW/GSM1036480_EB5K_01.txt"),
        header = TRUE
    )
    ensemble_ids <- counts$id
    df <- data.frame(row.names = ensemble_ids)

    # set columns with fpkms
    files <- list.files(
        path =
            "./benchmarkData/GSE42268/GSE42268_RAW",
        pattern = ".*.txt"
    )

    "
        gsm          cell_type1
        GSM1036483   G1
        ...
    "
    # not all cells have phases
    col_data_xml <- readRDS("./benchmarkData/ColDataXML.rds")
    gsms <- rownames(col_data_xml)
    for (file in files) {
        gsm <- substr(file, 1, 10)
        if (is.element(gsm, gsms)) {
            counts_data <- read.table(
                str_interp("./benchmarkData/GSE42268/GSE42268_RAW/${file}"),
                header = TRUE
            )
            counts_data <- subset(counts_data,
                select = -c(
                    id,
                    gene.symbol
                )
            )
            df[, gsm] <- counts_data
        }
    }

    # wilcox.test(#vector of all the values of the cells for each gene that are G2M, #same but for S/G1)
    # g2m <- c()
    # g1 <- c()
    # s <- c()
    gsmsG2M <- rownames(col_data_xml[col_data_xml$cell_type1 == "G2/M", ,
        drop = FALSE
    ])
    gsmsG1 <- rownames(col_data_xml[col_data_xml$cell_type1 == "G1", ,
        drop = FALSE
    ])
    gsmsS <- rownames(col_data_xml[col_data_xml$cell_type1 == "S", ,
        drop = FALSE
    ])

    countsG2M <- df[gsmsG2M]
    countsG1 <- df[gsmsG1]
    countsS <- df[gsmsS]
}

write_wilcox_output <- function(countsG2M, countsG1) {
    print("Writing to file...")
    file <- file("./output/wilcox.txt", "w")
    writeLines("id p_value", file)
    for (row in rownames(countsG2M)) {
        x <- countsG2M[row, ]
        y <- countsG1[row, ]
        # also display the mean and output to a file so that we can sort it by pvalue
        pval <- p.adjust(wilcox.test(as.numeric(unlist((x))), as.numeric(unlist(y)))$p.value, method = "fdr")
        line <- paste(row, pval, sep = " ")
        writeLines(line, file)
    }
    close(file)
    print("Writing done.")
}

read_wilcox_output <- function() {
    df <- read.table(
        "./output/wilcox.txt",
        header = TRUE
    )
    # how many are below 0.05
    print(paste0(nrow(df[df$p_value < 0.05, ]), " genes have a p-value below 0.05 out of ", nrow(df)))
    print(paste0(mean(df$p_value), " is the mean p-value"))
}

wilcox_fast_gse_42268()

read_wilcox_output()
