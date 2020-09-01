# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule busparse:
    input:
        gtf = "results/genomepy/{genome}/{genome}.annotation.gtf",
        fas = "results/genomepy/{genome}/{genome}.fa"
    output:
        ext = ["results/busparse/{genome}/cDNA_introns.fa",
               "results/busparse/{genome}/cDNA_tx_to_capture.txt",
               "results/busparse/{genome}/introns_tx_to_capture.txt",
               "results/busparse/{genome}/tr2g.tsv"]
    params:
        dir = "results/busparse/{genome}"
    log:
        "results/busparse/{genome}/get_velocity_files.log"
    message:
        "[busparse] Get files required for RNA velocity with bustools"
    conda:
        "../envs/busparse.yaml"
    script:
        "../scripts/busparse.R 2> {log}"