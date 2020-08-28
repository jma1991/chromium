#!/usr/bin/env Rscript

main <- function(input, output) {

    library(alevinQC)

    checkAlevinInputFiles(baseDir = input$dir)
    
    alevinQCReport(
        baseDir = input$dir,
        sampleId = wildcards$sample,
        outputFile = basename(output$pdf),
        outputFormat = "pdf_document",
        outputDir = dirname(output$pdf),
        forceOverwrite = TRUE
    )

}

main(snakemake@input, snakemake@output, snakemake@wildcards)