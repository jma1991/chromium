# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule gffread_tx2gene:
    input:
        "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        "results/gffread/{genome}/{genome}.tx2gene.tsv"
    message:
        "[gffread] Output a tx2gene annotation table"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input} --table transcript_id,gene_id | sort -u > {output}"

rule gffread_id2name:
    input:
        "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        "results/gffread/{genome}/{genome}.id2name.tsv"
    message:
        "[gffread] Output a id2name annotation table"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input} --table gene_id,gene_name | sort -u > {output}"

rule gffread_mrna:
	input:
		"results/genomepy/{genome}/{genome}.annotation.gtf"
	output:
		"results/gffread/{genome}/{genome}.mrna.txt"
	shell:
		"gffread {input} --table @chr,gene_id | grep MT | cut -f 2 | sort -u > {output}"

rule gffread_rrna:
	input:
		"results/genomepy/{genome}/{genome}.annotation.gtf"
	output:
		"results/gffread/{genome}/{genome}.rrna.txt"
	shell:
		"gffread {input} --table gene_biotype,gene_id | grep rRNA | cut -f 2 | sort -u > {output}"