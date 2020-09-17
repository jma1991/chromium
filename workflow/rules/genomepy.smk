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
    log:
        "results/genomepy/{genome}/{genome}.log"
    message:
        "[genomepy] Download reference genome and annotation: {wildcards.genome}"
    conda:
        "../envs/genomepy.yaml"
    shell:
        "genomepy install -g results/genomepy -a {wildcards.genome} 2> {log}"

rule genomepy_gunzip:
    input:
        "results/genomepy/{genome}/{genome}.annotation.{ext}.gz"
    output:
        "results/genomepy/{genome}/{genome}.annotation.{ext}"
    log:
        "results/genomepy/{genome}/{genome}.annotation.{ext}.log"
    message:
        "[genomepy] Extract reference annotation: {input}"
    conda:
        "../envs/genomepy.yaml"
    shell:
        "gunzip {input} 2> {log}"

rule genomepy_gentrome:
    input:
        "results/eisar/{genome}/{genome}.fa",
        "results/genomepy/{genome}/{genome}.fa"
    output:
        "results/genomepy/{genome}/{genome}.gentrome.fa"
    shell:
        "cat {input} > {output}"

rule genomepy_decoys:
	input:
		"results/genomepy/{genome}/{genome}.fa.fai"
	output:
		"results/genomepy/{genome}/{genome}.decoys.txt"
	shell:
		"cut -f 1 {input} > {output}"
