# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule singlecellexperiment_kallisto:
    input:
        rds = "results/busparse/read_velocity_output/{sample}.rds",
        tsv = expand("results/gffread/{genome}/{genome}.id2name.tsv", genome = config["genome"])
    output:
        rds = "results/singlecellexperiment/kallisto/{sample}.rds"
    params:
        out = "results/bustools/count/{sample}"
    log:
        out = "results/singlecellexperiment/kallisto/{sample}.out",
        err = "results/singlecellexperiment/kallisto/{sample}.err"
    message:
        "[singlecellexperiment] Create a SingleCellExperiment object from Kallisto output directory: {params.out}"
    conda:
        "../envs/singlecellexperiment.yaml"
    script:
        "../scripts/singlecellexperiment_kallisto.R"

rule singlecellexperiment_salmon:
    input:
        dir = "results/salmon/alevin/{sample}",
        tsv = [expand("results/eisar/{genome}/{genome}.features.tsv", genome = config["genome"]),
               expand("results/gffread/{genome}/{genome}.id2name.tsv", genome = config["genome"])]
    output:
        rds = "results/singlecellexperiment/salmon/{sample}.rds"
    log:
        out = "results/singlecellexperiment/salmon/{sample}.out",
        err = "results/singlecellexperiment/salmon/{sample}.err"
    message:
        "[singlecellexperiment] Create a SingleCellExperiment object from Salmon output directory: {input.dir}"
    conda:
        "../envs/singlecellexperiment.yaml"
    script:
        "../scripts/singlecellexperiment_salmon.R"

rule singlecellexperiment_star:
    input:
        dir = "results/star/align/{sample}"
    output:
        rds = "results/singlecellexperiment/star/{sample}.rds"
    log:
        out = "results/singlecellexperiment/star/{sample}.out",
        err = "results/singlecellexperiment/star/{sample}.err"
    message:
        "[singlecellexperiment] Create a SingleCellExperiment object from STAR output directory: {input.dir}"
    conda:
        "../envs/singlecellexperiment.yaml"
    script:
        "../scripts/singlecellexperiment_star.R"
