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

    library(eisaR)

    grl <- readRDS(input$rds)

    getTx2Gene(grl, filepath = output$tsv)

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)