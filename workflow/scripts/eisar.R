#!/usr/bin/env Rscript

main <- function(input, output) {

    library(Biostrings)

    library(eisaR)

    library(GenomicFeatures)

    library(Rsamtools)

    rng <- getFeatureRanges(
        gtf = input$gtf,
        featureType = c("spliced", "intron"),
        intronType = "separate",
        flankLength = 91,
        joinOverlappingIntrons = FALSE,
        verbose = TRUE
    )

    exportToGtf(rng, filepath = output$gtf)
    
    seq <- extractTranscriptSeqs(x = FaFile(input$fas), transcripts = rng)

    writeXStringSet(seq, filepath = output$fas)

    dat <- getTx2Gene(rng, filepath = output$tsv)

}

main(snakemake@input, snakemake@output)