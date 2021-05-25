# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule salmon_gentrome:
    input:
        "results/eisar/{genome}/{genome}.fa",
        "results/genomepy/{genome}/{genome}.fa"
    output:
        "results/salmon/genome/{genome}/gentrome.fa"
    log:
        "results/salmon/genome/{genome}/gentrome.log"
    conda:
        "../envs/salmon.yaml"
    shell:
        "cat {input} > {output} 2> {log}"

rule salmon_decoys:
    input:
        "results/genomepy/{genome}/{genome}.fa.fai"
    output:
        "results/salmon/genome/{genome}/decoys.txt"
    log:
        "results/salmon/genome/{genome}/decoys.log"
    conda:
        "../envs/salmon.yaml"
    shell:
        "cut -f 1 {input} > {output} 2> {log}"

rule salmon_index:
    input:
        transcripts = "results/salmon/genome/{genome}/gentrome.fa",
        decoys = "results/salmon/genome/{genome}/decoys.txt"
    output:
        index = directory("results/salmon/index/{genome}")
    log:
        out = "results/salmon/index/{genome}/log.out",
        err = "results/salmon/index/{genome}/log.err"
    message:
        "[salmon] Create a Salmon index"
    threads:
        16
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon index"
        " --transcripts {input.transcripts}"
        " --index {output.index}"
        " --threads {threads}"
        " --decoys {input.decoys}"
        " --keepDuplicates"
        " 1> {log.out}"
        " 2> {log.err}"

rule salmon_alevin:
    input:
        index = expand("results/salmon/index/{genome}", genome = config["genome"]),
        mates1 = lambda wildcards: units.loc[wildcards.sample, "read1"],
        mates2 = lambda wildcards: units.loc[wildcards.sample, "read2"],
        tgMap = expand("results/eisar/{genome}/{genome}.tx2gene.tsv", genome = config["genome"]),
        mrna = expand("results/gffread/{genome}/{genome}.mrna.txt", genome = config["genome"]),
        rrna = expand("results/gffread/{genome}/{genome}.rrna.txt", genome = config["genome"])
    output:
        output = directory("results/salmon/alevin/{sample}")
    params:
        arg = salmon_alevin_params
    log:
        out = "results/salmon/alevin/{sample}/log.out",
        err = "results/salmon/alevin/{sample}/log.err"
    message:
        "[salmon] Estimate gene abundances from scRNA-seq data"
    threads:
        16
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon alevin"
        " --libType ISR"
        " --index {input.index}"
        " --mates1 {input.mates1}"
        " --mates2 {input.mates2}"
        " {params.arg[protocol]}"
        " --output {output.output}"
        " --threads {threads}"
        " --tgMap {input.tgMap}"
        " --mrna {input.mrna}"
        " --rrna {input.rrna}"
        " 1> {log.out}"
        " 2> {log.err}"