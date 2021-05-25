# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, output, log, wildcards) {

     # Log function

    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    # Script function

    library(SingleCellExperiment)

    mat <- readRDS(input$rds)

    ann <- read.delim(input$tsv, header = FALSE, col.names = c("id", "name"))

    dim <- list(i = rownames(mat$spliced), j = colnames(mat$spliced))

    ind <- match(dim$i, ann$id)

    ann <- ann[ind, ]

    sce <- SingleCellExperiment(
        
        assays = list(
            counts = mat$spliced,
            spliced = mat$spliced,
            unspliced = mat$unspliced
        ),
        
        rowData = DataFrame(
            ID = ann$id,
            Symbol = ann$name
        ),
        
        colData = DataFrame(
            Sample = wildcards$sample,
            Barcode = dim$j
        )
    )

    saveRDS(sce, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@log, snakemake@wildcards)