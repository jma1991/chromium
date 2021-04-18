# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, output, log) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    library(GenomicFeatures)

    library(Rsamtools)

    grl <- readRDS(input$rds)

    con <- FaFile(input$fas)

    seq <- extractTranscriptSeqs(con, grl)

    writeXStringSet(seq, output$fas)

}

main(snakemake@input, snakemake@output, snakemake@log)