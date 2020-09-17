# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule salmon_index:
    input:
        fas = "results/genomepy/{genome}/{genome}.gentrome.fa",
        txt = "results/genomepy/{genome}/{genome}.decoys.txt"
    output:
        idx = directory("results/salmon/{genome}")
    log:
        out = "results/salmon/{genome}/log.out",
        err = "results/salmon/{genome}/log.err"
    message:
        "[salmon] Create a Salmon index"
    threads:
        16
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon index -t {input.fas} -i {output.idx} -p {threads} -d {input.txt} --keepDuplicates 1> {log.out} 2> {log.err}"

rule salmon_alevin:
    input:
        idx = expand("results/salmon/{genome}", genome = config["genome"]),
        fq1 = lambda wildcards: expand("results/cutadapt/{{sample}}_{unit}_1.fastq.gz", unit = units.loc[wildcards.sample, "unit"]),
        fq2 = lambda wildcards: expand("results/cutadapt/{{sample}}_{unit}_2.fastq.gz", unit = units.loc[wildcards.sample, "unit"]),
        tsv = expand("results/eisar/{genome}/{genome}.tx2gene.tsv", genome = config["genome"]),
        txt = expand(["results/gffread/{genome}/{genome}.mrna.txt", "results/gffread/{genome}/{genome}.rrna.txt"], genome = config["genome"])
    output:
        mat = "results/salmon/{sample}/alevin/quants_mat.gz"
    log:
        out = "results/salmon/{sample}/log.out",
        err = "results/salmon/{sample}/log.err"
    params:
        out = "results/salmon/{sample}"
    message:
        "[salmon] Estimate gene abundances from scRNA-seq data"
    threads:
        16
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon alevin -l ISR -i {input.idx} -1 {input.fq1} -2 {input.fq2} -o {params.out} -p {threads} --tgMap {input.tsv} --chromiumV3 --mrna {input.txt[0]} --rrna {input.txt[1]} --dumpFeatures 1> {log.out} 2> {log.err}"
