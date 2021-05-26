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

    len <- switch(params$chemistry, "10xv1" = 98, "10xv2" = 98, "10xv3" = 91)

    grl <- getFeatureRanges(
        gtf = input$gtf,
        featureType = c("spliced", "intron"),
        intronType = "separate",
        flankLength = len,
        joinOverlappingIntrons = FALSE,
        verbose = TRUE
    )

    saveRDS(grl, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log)