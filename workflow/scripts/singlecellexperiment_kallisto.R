#!/usr/bin/env Rscript

main <- function(input, output, wildcards) {

    # Load Bioconductor packages
    
    library(BUSpaRse)
    
    library(SingleCellExperiment)

    #

    spliced_mtx <- sub('\\.mtx$', '', input$mtx[1])

    unspliced_mtx <- sub('\\.mtx$', '', input$mtx[2])

    out <- read_velocity_output(
        spliced_dir = dirname(spliced_mtx),
        spliced_name = basename(spliced_mtx),
        unspliced_dir = dirname(unspliced_mtx),
        unspliced_name = basename(unspliced_mtx)
    )
    
    ann <- read.delim(input$tsv, header = FALSE, col.names = c("id", "name"))
    
    use <- intersect(colnames(out$spliced), colnames(out$unspliced))
    
    out$spliced <- out$spliced[ann$id, use]
    
    out$unspliced <- out$unspliced[ann$id, use]
    
    sce <- SingleCellExperiment(
        assays = list(counts = out$spliced, spliced = out$spliced, unspliced = out$unspliced),
        rowData = DataFrame(ID = ann$id, Symbol = ann$name),
        colData = DataFrame(Sample = wildcards$sample, Barcode = use)
    )

    saveRDS(sce, file = output$rds)

}

main(snakemake@input, snakemake@output, snakemake@wildcards)