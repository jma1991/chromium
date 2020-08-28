#!/usr/bin/env Rscript

main <- function(input, output) {
    
    library(Matrix)
    
    
    
    readLines(con = "")

    genes <- readr::read_delim("Velocyto/raw/features.tsv", col_names = FALSE, delim = "\t")
    
    cells <- readr::read_delim("Velocyto/raw/barcodes.tsv", col_names = FALSE, delim = "\t")
    
    n1 <- nrow(genes)
    
    n2 <- nrow(cells)
    
    ## Read matrix. col 1 = gene, col 2 = cell, col 3 = spliced, col 4 = unspliced, col 5 = ambiguous
    matr <- readr::read_delim("Velocyto/raw/matrix.mtx", skip = 3, delim = " ", col_names = FALSE)
    
    scounts <- sparseMatrix(i = m$X1, j = m$X2, x = m$X3, dims = c(n1, n2))
    
    ucounts <- sparseMatrix(i = m$X1, j = m$X2, x = m$X4, dims = c(n1, n2))
    
    rownames(scounts) <- genes$X2
    
    rownames(ucounts) <- genes$X2
    
    colnames(scounts) <- cells$X1
        
    colnames(ucounts) <- cells$X1
    
    sce <- SingleCellExperiment(
        assays = list(counts = scounts,
                      spliced = scounts,
                      unspliced = ucounts)
    )
    
    
}

main(snakemake@input, snakemake@output)