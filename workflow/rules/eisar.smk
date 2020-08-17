# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule eisar:
    input:
        gtf = "results/genomepy/{genome}/{genome}.annotation.gtf",
        fas = "results/genomepy/{genome}/{genome}.fa"
    output:
        gtf = "results/eisar/{genome}/{genome}.annotation.expanded.gtf",
        fas = "results/eisar/{genome}/{genome}.annotation.expanded.fa",
        tsv = "results/eisar/{genome}/{genome}.annotation.expanded.tx2gene.tsv"
    script:
        "../scripts/eisar.R"
