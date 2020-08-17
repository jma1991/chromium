#!/usr/bin/env Rscript

main <- function(input, output) {
    
    library(Matrix)

    library(SingleCellExperiment)
    
    mat <- t(readMM(input$mat))

    row <- DataFrame(
        ID = readLines(input$row),
        Symbol = NULL
    )

    rownames(mat) <- row$ID
    
    col <- DataFrame(
        Sample = NULL,
        Barcode = readLines(input$col),
    )

    colnames(mat) <- col$Barcode
        
    sce <- SingleCellExperiment(
        assays = list(counts = t(mat)),
        rowData = row,
        colData = col
    )

    saveRDS(sce, file = output$sce)

}

main(snakemake@input, snakemake@output)