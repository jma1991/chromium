# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

def kallisto_technology(wildcards):
    arg = ["10xv1", "10xv2", "10xv3", "CELSeq", "CELSeq2", "DropSeq", "inDropsv1", "inDropsv2", "inDropsv3", "SCRBSeq", "SureCell"]

rule kallisto_index:
    input:
        fas = "results/genomepy/{genome}/{genome}.fa"
    output:
        idx = "results/genomepy/{genome}/{genome}/index/kallisto/{genome}.idx"
    message:
        "[kallisto] Build kallisto index"
    conda:
        "../rules/kallisto.yaml"
    shell:
        "kallisto index -i {output.idx} -k 31 {input.fas}"

rule kallisto_bus:
    input:
        idx = "results/genomepy/{genome}/{genome}/index/kallisto/{genome}.idx"
    output:
        dir = directory("results/kallisto/{sample}")
    params:
        arg = ""
    threads:
        16
    message:
        "[kallisto] Generate BUS file for single-cell data"
    conda:
        "../rules/kallisto.yaml"
    shell:
        "kallisto bus -i {input.idx} -o {output.dir} -x {params.arg} -t {threads} {input.fas}"
