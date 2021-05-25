# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule gffread_tx2gene:
    input:
        "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        "results/gffread/{genome}/{genome}.tx2gene.tsv"
    log:
        "results/gffread/{genome}/{genome}.tx2gene.log"
    message:
        "[gffread] Create transcript_id to gene_id annotation table"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input} --table transcript_id,gene_id | sort -u > {output} 2> {log}"

rule gffread_id2name:
    input:
        "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        "results/gffread/{genome}/{genome}.id2name.tsv"
    log:
        "results/gffread/{genome}/{genome}.id2name.log"
    message:
        "[gffread] Create gene_id to gene_name annotation table"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input} --table gene_id,gene_name | sort -u > {output} 2> {log}"

rule gffread_mrna:
    input:
        "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        "results/gffread/{genome}/{genome}.mrna.txt"
    log:
        "results/gffread/{genome}/{genome}.mrna.log"
    message:
        "[gffread] Create mRNA annotation table"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input} --table @chr,gene_id | awk '$1 == MT' | cut -f 2 | sort -u > {output} 2> {log}"

rule gffread_rrna:
    input:
        "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        "results/gffread/{genome}/{genome}.rrna.txt"
    log:
        "results/gffread/{genome}/{genome}.rrna.log"
    message:
        "[gffread] Create rRNA annotation table"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input} --table gene_biotype,gene_id | awk '$1 == rRNA' | cut -f 2 | sort -u > {output} 2> {log}"
