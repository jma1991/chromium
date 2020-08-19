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

    # Write spliced and unspliced gene identifiers to disk
    write.table(x = metadata(rng)$corrgene, file = output$tsv[1], row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    # Write transcript and intron identifiers to disk
    getTx2Gene(rng, filepath = output$tsv[2])

}

main(snakemake@input, snakemake@output)