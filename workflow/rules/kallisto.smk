# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule kallisto_index:
    input:
        fas = "results/busparse/get_velocity_files/{genome}/cDNA_introns.fa"
    output:
        idx = "results/kallisto/index/{genome}/{genome}.idx"
    log:
        log = "results/kallisto/index/{genome}/{genome}.log"
    message:
        "[kallisto] Build kallisto index"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index"
        " --index {output.idx}"
        " {input.fas}"
        " 2> {log}"

rule kallisto_bus:
    input:
        idx = expand("results/kallisto/index/{genome}/{genome}.idx", genome = config["genome"]),
        fq1 = lambda wildcards: units.loc[wildcards.sample, "read1"],
        fq2 = lambda wildcards: units.loc[wildcards.sample, "read2"]
    output:
        ext = ["results/kallisto/bus/{sample}/matrix.ec",
               "results/kallisto/bus/{sample}/output.bus",
               "results/kallisto/bus/{sample}/run_info.json",
               "results/kallisto/bus/{sample}/transcripts.txt"]
    params:
        out = lambda wildcards, output: os.path.dirname(output[0]),
        ver = config["chemistry"],
        fqz = lambda wildcards, input: [j for i in zip(input.fq1, input.fq2) for j in i]
    log:
        log = "results/kallisto/bus/{sample}/kallisto_bus.log"
    message:
        "[kallisto] Generate BUS file for single-cell data"
    threads:
        16
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto bus"
        " --index {input.idx}"
        " --output-dir {params.out}"
        " --technology {params.ver}"
        " --threads {threads}"
        " {params.fqz}"
        " 2> {log}"