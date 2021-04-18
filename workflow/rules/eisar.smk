# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule eisar_ranges:
    input:
        gtf = "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        rds = "results/eisar/{genome}/{genome}.ranges.rds"
    params:
        chemistry = config["chemistry"]
    log:
        out = "results/eisar/{genome}/{genome}.ranges.out",
        err = "results/eisar/{genome}/{genome}.ranges.err"
    message:
        "[eisar] Generate a GRangesList object with feature ranges"
    conda:
        "../envs/eisar.yaml"
    script:
        "../scripts/eisar_ranges.R"

rule eisar_fa:
    input:
        rds = "results/eisar/{genome}/{genome}.ranges.rds",
        fas = "results/genomepy/{genome}/{genome}.fa"
    output:
        fas = "results/eisar/{genome}/{genome}.fa"
    log:
        out = "results/eisar/{genome}/{genome}.out",
        err = "results/eisar/{genome}/{genome}.err"
    message:
        "[eisar] Extract transcript (or CDS) sequences from chromosome sequences"
    conda:
        "../envs/eisar.yaml"
    script:
        "../scripts/eisar_fa.R"

rule eisar_gtf:
    input:
        rds = "results/eisar/{genome}/{genome}.ranges.rds"
    output:
        gtf = "results/eisar/{genome}/{genome}.annotation.gtf"
    log:
        out = "results/eisar/{genome}/{genome}.annotation.out",
        err = "results/eisar/{genome}/{genome}.annotation.err"
    message:
        "[eisar] Export GRangesList to GTF"
    conda:
        "../envs/eisar.yaml"
    script:
        "../scripts/eisar_gtf.R"

rule eisar_tsv:
    input:
        rds = "results/eisar/{genome}/{genome}.ranges.rds"
    output:
        tsv = "results/eisar/{genome}/{genome}.features.tsv"
    log:
        out = "results/eisar/{genome}/{genome}.features.out",
        err = "results/eisar/{genome}/{genome}.features.err"
    message:
        "[eisar] Generate a spliced-to-intron mapping from a GRangesList"
    conda:
        "../envs/eisar.yaml"
    script:
        "../scripts/eisar_tsv.R"

rule eisar_tx2gene:
    input:
        rds = "results/eisar/{genome}/{genome}.ranges.rds"
    output:
        tsv = "results/eisar/{genome}/{genome}.tx2gene.tsv"
    log:
        out = "results/eisar/{genome}/{genome}.tx2gene.out",
        err = "results/eisar/{genome}/{genome}.tx2gene.err"
    message:
        "[eisar] Generate a transcript-to-gene mapping from a GRangesList"
    conda:
        "../envs/eisar.yaml"
    script:
        "../scripts/eisar_tx2gene.R"
