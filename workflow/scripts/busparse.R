#!/usr/bin/env Rscript

main <- function(input, params) {

    library(Biostrings)

    library(BUSpaRse)
    
    get_velocity_files(
        X = input$gtf,
        L = 91,
        Genome = readDNAStringSet(input$fas),
        out_path = params$dir,
        style = "Ensembl"
    )

}

main(snakemake@input, snakemake@params)