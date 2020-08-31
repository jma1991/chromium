#!/usr/bin/env Rscript

main <- function(input, output, params, wildcards) {


    # Read barcodes

    barcodes <- file.path(params$out, "barcodes.tsv")

    names <- c("barcode")

    barcodes <- data.table::fread(barcodes, header = FALSE, col.names = names)


    # Read features

    features <- file.path(params$out, "features.tsv")

    names  <- c("id", "name", "type")

    features <- data.table::fread(features, header = FALSE, col.names = names)


    # Read matrix

    matrix <- file.path(params$out, "matrix.mtx")

    names <- c("feature", "barcode", "spliced", "unspliced", "ambiguous")

    matrix <- data.table::fread(matrix, header = FALSE, skip = 3, col.names = names)


    # Calculate dimensions

    shape <- c(nrow(features), nrow(barcodes))

    names <- list(features$id, barcodes$barcode)


    # Create object

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
            counts = Matrix::sparseMatrix(i = matrix$feature, j = matrix$barcode, x = matrix$spliced, dims = shape, dimnames = names),
            spliced = Matrix::sparseMatrix(i = matrix$feature, j = matrix$barcode, x = matrix$spliced, dims = shape, dimnames = names),
            unspliced = Matrix::sparseMatrix(i = matrix$feature, j = matrix$barcode, x = matrix$unspliced, dims = shape, dimnames = names)
        ),
        rowData = S4Vectors::DataFrame(
            ID = features$id,
            Symbol = features$name
        ),
        colData = S4Vectors::DataFrame(
            Sample = wildcards$sample,
            Barcode = barcodes$barcode
        )
    )


    # Save object

    saveRDS(sce, file = output$rds)
    
}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@wildcards)
