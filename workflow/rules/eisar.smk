# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule eisar:
    input:
        fas = "results/genomepy/{genome}/{genome}.fa",
        gtf = "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        fas = "results/eisar/{genome}/{genome}.fa",
        gtf = "results/eisar/{genome}/{genome}.annotation.gtf",
        tsv = ["results/eisar/{genome}/{genome}.features.tsv", "results/eisar/{genome}/{genome}.tx2gene.tsv"]
    script:
        "../scripts/eisar.R"
