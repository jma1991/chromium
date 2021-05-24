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

    library(BUSpaRse)

    spliced_mtx <- sub('\\.mtx$', '', input$mtx[1])

    unspliced_mtx <- sub('\\.mtx$', '', input$mtx[2])

    out <- read_velocity_output(
        spliced_dir = dirname(spliced_mtx),
        unspliced_dir = dirname(unspliced_mtx),
        spliced_name = basename(spliced_mtx),
        unspliced_name = basename(unspliced_mtx)
    )

    row <- intersect(rownames(out$spliced), rownames(out$unspliced))

    col <- intersect(colnames(out$spliced), colnames(out$unspliced))

    out$spliced <- out$spliced[row, col]

    out$unspliced <- out$unspliced[row, col]

    saveRDS(out, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@log)