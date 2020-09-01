#!/usr/bin/env Rscript

main <- function(input, output, params, wildcards) {

    library(alevinQC)

    checkAlevinInputFiles(baseDir = params$dir)
    
    alevinQCReport(
        baseDir = params$dir,
        sampleId = wildcards$sample,
        outputFile = basename(output$html),
        outputFormat = "html_document",
        outputDir = dirname(output$html),
        forceOverwrite = TRUE
    )

}

main(snakemake@input, snakemake@output, snakemake@params, snakemake@wildcards)