# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule singlecellexperiment_kallisto:
    input:
        mtx = expand("results/bustools/{{sample}}/output.correct.sort.{feature}.mtx", feature = ["spliced", "unspliced"]),
        tsv = "results/gffread/GRCm38.p6/GRCm38.p6.id2name.tsv"
    output:
        rds = "results/singlecellexperiment/{sample}/kallisto.rds"
    params:
        out = "results/bustools/{sample}"
    message:
        "[singlecellexperiment] Create a SingleCellExperiment object from Kallisto output directory: {params.out}"
    conda:
        "../envs/singlecellexperiment.yaml"
    script:
        "../scripts/singlecellexperiment_kallisto.R"

rule singlecellexperiment_salmon:
    input:
        mat = "results/salmon/alevin/{sample}/alevin/quants_mat.gz",
        tsv = ["results/eisar/GRCm38.p6/GRCm38.p6.features.tsv", "results/gffread/GRCm38.p6/GRCm38.p6.id2name.tsv"]
    output:
        rds = "results/singlecellexperiment/{sample}/salmon.rds"
    params:
        out = "results/salmon/alevin/{sample}/alevin"
    message:
        "[singlecellexperiment] Create a SingleCellExperiment object from Salmon output directory: {params.out}"
    conda:
        "../envs/singlecellexperiment.yaml"
    script:
        "../scripts/singlecellexperiment_salmon.R"

rule singlecellexperiment_star:
    input:
        mtx = "results/star/solo/{sample}/Solo.out/Velocyto/raw/matrix.mtx"
    output:
        rds = "results/singlecellexperiment/{sample}/star.rds"
    params:
        out = "results/star/solo/{sample}/Solo.out/Velocyto/raw"
    message:
        "[singlecellexperiment] Create a SingleCellExperiment object from STAR output directory: {params.out}"
    conda:
        "../envs/singlecellexperiment.yaml"
    script:
        "../scripts/singlecellexperiment_star.R"