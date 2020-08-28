#!/usr/bin/env Rscript

main <- function(input, output, wildcards) {

    # Load Bioconductor packages

    library(tximport)

    library(tximeta)

    library(SummarizedExperiment)

    library(SingleCellExperiment)

    # Import salmon quantification

    txi <- tximport(input$mat, type = "alevin", dropInfReps = FALSE)
    
    rse <- SummarizedExperiment(assays = list(counts = round(txi$counts)))

    # Split by spliced and unspliced
    
    ann <- read.delim(input$tsv[1], header = TRUE, col.names = c("spliced", "unspliced"))
    
    rse <- splitSE(rse, ann, assayName = "counts")

    # Import and match annotation

    ann <- read.delim(input$tsv[2], header = FALSE, col.names = c("id", "name"))

    rse <- rse[ann$id, ]

    # Create expanded SCE object

    sce <- SingleCellExperiment(
        assays = list(
            counts = assay(rse, "spliced"),
            spliced = assay(rse, "spliced"),
            unspliced = assay(rse, "unspliced")
        ),
        rowData = DataFrame(
            ID = rownames(rse),
            Symbol = ann$name
        ),
        colData = DataFrame(
            Sample = wildcards$sample,
            Barcode = colnames(rse)
        )
    )

    # Save SCE object to disk

    saveRDS(sce, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@wildcards)