# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, output, params, log) {

    # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    library(tximeta)

    makeLinkedTxome(
        indexDir = input$indexDir,
        source = params$source,
        organism = params$organism,
        release = params$release,
        genome = params$genome,
        fasta = input$fasta,
        gtf = input$gtf,
        jsonFile = output$jsonFile
    )

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)