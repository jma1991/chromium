#!/usr/bin/env Rscript

main <- function(input, log, params) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

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
        gene_version = NULL,
        chrs_only = FALSE
    )

}

main(snakemake@input, snakemake@log, snakemake@params)
