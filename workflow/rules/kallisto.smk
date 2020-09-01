# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule kallisto_index:
    input:
        fas = "results/busparse/{genome}/cDNA_introns.fa"
    output:
        idx = "results/kallisto/{genome}.idx"
    log:
        log = "results/kallisto/{genome}.log"
    message:
        "[kallisto] Build kallisto index"
    conda:
        "../rules/kallisto.yaml"
    shell:
        "kallisto index -i {output.idx} {input.fas} 2> {log}"

rule kallisto_bus:
    input:
        idx = "results/kallisto/GRCm38.p6.idx",
        fq1 = lambda wildcards: pep.subsample_table.loc[pep.subsample_table['sample_name'] == wildcards.sample, "read1"],
        fq2 = lambda wildcards: pep.subsample_table.loc[pep.subsample_table['sample_name'] == wildcards.sample, "read2"]
    output:
        ext = ["results/kallisto/{sample}/matrix.ec",
               "results/kallisto/{sample}/output.bus",
               "results/kallisto/{sample}/run_info.json",
               "results/kallisto/{sample}/transcripts.txt"]
    log:
        log = "results/kallisto/{sample}/kallisto_bus.log"
    params:
        out = "results/kallisto/{sample}",
        fqz = lambda wildcards, input: [j for i in zip(input.fq1, input.fq2) for j in i]
    threads:
        16
    message:
        "[kallisto] Generate BUS file for single-cell data"
    conda:
        "../rules/kallisto.yaml"
    shell:
        "kallisto bus -i {input.idx} -o {params.out} -x 10xv3 -t {threads} {params.fqz} 2> {log}"
