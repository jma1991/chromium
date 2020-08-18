
genes <- readr::read_delim("Velocyto/raw/features.tsv", col_names = FALSE, delim = "\t")

cells <- readr::read_delim("Velocyto/raw/barcodes.tsv", col_names = FALSE, delim = "\t")

n_genes <- nrow(genes)

n_cells <- nrow(cells)

## Read matrix. col 1 = gene, col 2 = cell, col 3 = spliced, col 4 = unspliced, col 5 = ambiguous
matr <- readr::read_delim("Velocyto/raw/matrix.mtx", skip = 3, delim = " ", col_names = FALSE)

scounts <- Matrix::sparseMatrix(
    i = matr$X1,
    j = matr$X2,
    x = matr$X3,
    dims = c(n_genes, n_cells)
)

ucounts <- Matrix::sparseMatrix(
    i = matr$X1,
    j = matr$X2,
    x = matr$X4,
    dims = c(n_genes, n_cells)
)

rownames(scounts) <- rownames(ucounts) <- genes$X2

colnames(scounts) <- colnames(ucounts) <- cells$X1

sce <- SingleCellExperiment(
        assays = list(counts = scounts,
                      spliced = scounts,
                      unspliced = ucounts)
    )