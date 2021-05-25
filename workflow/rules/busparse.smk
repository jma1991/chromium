# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule busparse_get_velocity_files:
    input:
        gtf = "results/genomepy/{genome}/{genome}.annotation.gtf",
        fas = "results/genomepy/{genome}/{genome}.fa"
    output:
        ext = ["results/busparse/get_velocity_files/{genome}/cDNA_introns.fa",
               "results/busparse/get_velocity_files/{genome}/cDNA_tx_to_capture.txt",
               "results/busparse/get_velocity_files/{genome}/introns_tx_to_capture.txt",
               "results/busparse/get_velocity_files/{genome}/tr2g.tsv"]
    params:
        chemistry = config["chemistry"],
        style = config["source"],
        out_path = lambda wildcards, output: os.path.dirname(output[0])
    log:
        out = "results/busparse/get_velocity_files/{genome}/log.out",
        err = "results/busparse/get_velocity_files/{genome}/log.err"
    message:
        "[busparse] Get files required for RNA velocity with bustools"
    conda:
        "../envs/busparse.yaml"
    script:
        "../scripts/get_velocity_files.R"

rule busparse_read_velocity_output:
    input:
        mtx = ["results/bustools/count/{sample}/spliced.mtx",
               "results/bustools/count/{sample}/unspliced.mtx"]
    output:
        rds = "results/busparse/read_velocity_output/{sample}.rds"
    log:
        out = "results/busparse/read_velocity_output/{sample}.out",
        err = "results/busparse/read_velocity_output/{sample}.err"
    message:
        "[busparse] Read intronic and exonic matrices into R"
    conda:
        "../envs/busparse.yaml"
    script:
        "../scripts/read_velocity_output.R"