#!/usr/bin/env Rscript

main <- function(input, output) {

    library(tximeta)
    
    makeLinkedTxome(
        indexDir = input$idx,
        source = "",
        genome = "",
        organism = "",
        release = "",
        fasta = input$fas,
        gtf = input$gtf,
        write = TRUE,
        jsonFile = output$jsn
    )

}

main(snakemake@input, snakemake@output)