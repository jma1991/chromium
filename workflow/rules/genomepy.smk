# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule genomepy_install:
    output:
        "results/genomepy/{genome}.annotation.bed.gz",
        "results/genomepy/{genome}.annotation.gtf.gz",
        "results/genomepy/{genome}.fa",
        "results/genomepy/{genome}.fa.fai",
        "results/genomepy/{genome}.fa.sizes",
        "results/genomepy/{genome}.gaps.bed",
        "results/genomepy/README.txt"
    params:
        "results/genomepy"
    log:
        "results/genomepy/{genome}.log"
    message:
        "[genomepy] Download reference genome and annotation: {wildcards.genome}"
    conda:
        "../envs/genomepy.yaml"
    shell:
        "genomepy install -g {params} -a {wildcards.genome} 2> {log}"

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
        "gunzip -c {input} > {output}"
