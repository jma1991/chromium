# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, output, log, wildcards) {


    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")


    # Load Bioconductor packages

    library(tximport)

    library(tximeta)

    library(SummarizedExperiment)

    library(SingleCellExperiment)


    # Import salmon quantification

    con <- file.path(input$dir, "alevin", "quants_mat.gz")

    txi <- tximport(con, type = "alevin", dropInfReps = FALSE)
    
    rse <- SummarizedExperiment(assays = list(counts = round(txi$counts)))


    # Split by spliced and unspliced
    
    dat <- read.delim(input$tsv[1], header = TRUE, col.names = c("spliced", "unspliced"))
    
    rse <- splitSE(rse, dat, assayName = "counts")


    # Import and match annotation

    ann <- read.delim(input$tsv[2], header = FALSE, col.names = c("id", "name"))

    ann <- ann[match(rownames(rse), ann$id), ]


    # Create SCE object

    sce <- SingleCellExperiment(
        assays = list(
            counts = assay(rse, "spliced"),
            spliced = assay(rse, "spliced"),
            unspliced = assay(rse, "unspliced")
        ),
        rowData = DataFrame(
            ID = ann$id,
            Symbol = ann$name
        ),
        colData = DataFrame(
            Sample = wildcards$sample,
            Barcode = colnames(rse)
        )
    )


    # Save SCE object

    saveRDS(sce, file = output$rds)


}

main(snakemake@input, snakemake@output, snakemake@log, snakemake@wildcards)
