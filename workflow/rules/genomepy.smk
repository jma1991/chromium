# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule genomepy_install:
    output:
        "results/genomepy/{genome}/{genome}.annotation.bed.gz",
        "results/genomepy/{genome}/{genome}.annotation.gtf.gz",
        "results/genomepy/{genome}/{genome}.fa",
        "results/genomepy/{genome}/{genome}.fa.fai",
        "results/genomepy/{genome}/{genome}.fa.sizes",
        "results/genomepy/{genome}/{genome}.gaps.bed",
        "results/genomepy/{genome}/README.txt"
    message:
        "[genomepy] Download reference genome and annotation: {wildcards.genome}"
    conda:
        "../envs/genomepy.yaml"
    shell:
        "genomepy install -g results/genomepy -a {wildcards.genome} Ensembl"

rule genomepy_gunzip:
    input:
        "results/genomepy/{genome}/{genome}.annotation.{ext}.gz"
    output:
        "results/genomepy/{genome}/{genome}.annotation.{ext}"
    message:
        "[genomepy] Extract reference annotation: {input}"
    conda:
        "../envs/genomepy.yaml"
    shell:
        "gunzip {input}"