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

    library(GenomicRanges)

    grl <- readRDS(input$rds)

    write.table(
        x = metadata(grl)$corrgene,
        file = output$tsv,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)