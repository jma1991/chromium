# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

main <- function(input, output, params, log, wildcards) {


    # Setup logging
    
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")


    # Define directory
    
    dirname <- input$dir

    dirname <- file.path(dirname, "Solo.out", "Velocyto", "raw")


    # Read features

    names  <- c("id", "name", "type")

    features <- file.path(dirname, "features.tsv")

    features <- data.table::fread(features, header = FALSE, col.names = names)


    # Read barcodes

    names <- c("id")

    barcodes <- file.path(dirname, "barcodes.tsv")

    barcodes <- data.table::fread(barcodes, header = FALSE, col.names = names)

    # Read spliced matrix

    names <- c("feature", "barcode", "count")

    spliced <- file.path(dirname, "spliced.mtx")

    spliced <- data.table::fread(spliced, header = FALSE, skip = 2, col.names = names)


    # Read unspliced matrix

    names <- c("feature", "barcode", "count")

    unspliced <- file.path(dirname, "unspliced.mtx")

    unspliced <- data.table::fread(unspliced, header = FALSE, skip = 2, col.names = names)


    # Set dimensions

    shape <- c(nrow(features), nrow(barcodes))

    names <- list(features$id, barcodes$id)


    # Create object

    sce <- SingleCellExperiment::SingleCellExperiment(
        
        assays = list(
            
            counts = Matrix::sparseMatrix(
                i = spliced$feature,
                j = spliced$barcode,
                x = spliced$count,
                dims = shape,
                dimnames = names
            ),
            
            spliced = Matrix::sparseMatrix(
                i = spliced$feature,
                j = spliced$barcode,
                x = spliced$count,
                dims = shape,
                dimnames = names
            ),
            
            unspliced = Matrix::sparseMatrix(
                i = unspliced$feature,
                j = unspliced$barcode,
                x = unspliced$count,
                dims = shape,
                dimnames = names
            )
        ),
        
        rowData = S4Vectors::DataFrame(
            ID = features$id,
            Symbol = features$name
        ),
        
        colData = S4Vectors::DataFrame(
            Sample = wildcards$sample,
            Barcode = barcodes$id
        )
    
    )

    # Save object

    saveRDS(sce, file = output$rds)


}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@log, snakemake@wildcards)
