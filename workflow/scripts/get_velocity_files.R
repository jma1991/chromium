# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, params, log) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    library(Biostrings)

    library(BUSpaRse)

    len <- switch(params$chemistry, "10xv1" = 98, "10xv2" = 98, "10xv3" = 91)

    dna <- readDNAStringSet(input$fas)

    names(dna) <- sapply(strsplit(names(dna), " "), .subset, 1)
    
    get_velocity_files(
        X = input$gtf,
        L = len,
        Genome = dna,
        out_path = params$out_path,
        style = params$style,
        transcript_version = NULL,
        gene_version = NULL
    )

}

main(snakemake@input, snakemake@params, snakemake@log)
