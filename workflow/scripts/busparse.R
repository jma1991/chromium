#!/usr/bin/env Rscript

main <- function(input, params) {

    library(Biostrings)

    library(BUSpaRse)

    dna <- readDNAStringSet(input$fas)

    names(dna) <- sapply(strsplit(names(dna), " "), .subset, 1)
    
    get_velocity_files(
        X = input$gtf,
        L = 91,
        Genome = dna,
        out_path = params$dir,
        style = "Ensembl",
        transcript_version = NULL,
        gene_version = NULL
    )

}

main(snakemake@input, snakemake@params)