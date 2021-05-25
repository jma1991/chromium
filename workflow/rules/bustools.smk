# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule bustools_correct:
    input:
        txt = expand("resources/barcodes/{chemistry}.txt", chemistry = config["chemistry"]),
        bus = "results/kallisto/bus/{sample}/output.bus"
    output:
        bus = "results/bustools/correct/{sample}/output.bus"
    log:
        log = "results/bustools/correct/{sample}/output.log"
    message:
        "[bustools] Error correct BUS file"
    conda:
        "../envs/bustools.yaml"
    shell:
        "bustools correct -o {output.bus} -w {input.txt} {input.bus} 2> {log}"

rule bustools_sort:
    input:
        bus = "results/bustools/correct/{sample}/output.bus"
    output:
        bus = "results/bustools/sort/{sample}/output.bus"
    log:
        log = "results/bustools/sort/{sample}/output.log"
    message:
        "[bustools] Sort BUS file by barcode and UMI"
    conda:
        "../envs/bustools.yaml"
    shell:
        "bustools sort -t {threads} -o {output.bus} {input.bus} 2> {log}"

rule bustools_capture_spliced:
    input:
        txt = expand("results/busparse/get_velocity_files/{genome}/introns_tx_to_capture.txt", genome = config["genome"]),
        mat = "results/kallisto/bus/{sample}/matrix.ec",
        txi = "results/kallisto/bus/{sample}/transcripts.txt",
        bus = "results/bustools/sort/{sample}/output.bus"
    output:
        bus = "results/bustools/capture/{sample}/spliced.bus"
    log:
        log = "results/bustools/capture/{sample}/spliced.log"
    message:
        "[bustools] Capture spliced reads from BUS file"
    conda:
        "../envs/bustools.yaml"
    shell:
        "bustools capture -s -x -o {output.bus} -c {input.txt} -e {input.mat} -t {input.txi} {input.bus} 2> {log}"

rule bustools_capture_unspliced:
    input:
        txt = expand("results/busparse/get_velocity_files/{genome}/cDNA_tx_to_capture.txt", genome = config["genome"]),
        mat = "results/kallisto/bus/{sample}/matrix.ec",
        txi = "results/kallisto/bus/{sample}/transcripts.txt",
        bus = "results/bustools/sort/{sample}/output.bus"
    output:
        bus = "results/bustools/capture/{sample}/unspliced.bus"
    log:
        log = "results/bustools/capture/{sample}/unspliced.log"
    message:
        "[bustools] Capture unspliced reads from BUS file"
    conda:
        "../envs/bustools.yaml"
    shell:
        "bustools capture -s -x -o {output.bus} -c {input.txt} -e {input.mat} -t {input.txi} {input.bus} 2> {log}"

rule bustools_count_spliced:
    input:
        tsv = expand("results/busparse/get_velocity_files/{genome}/tr2g.tsv", genome = config["genome"]),
        mat = "results/kallisto/bus/{sample}/matrix.ec",
        txt = "results/kallisto/bus/{sample}/transcripts.txt",
        bus = "results/bustools/capture/{sample}/spliced.bus"
    output:
        ext = multiext("results/bustools/count/{sample}/spliced", ".barcodes.txt", ".genes.txt", ".mtx")
    log:
        log = "results/bustools/count/{sample}/spliced.log"
    params:
        out = lambda wildcards, output: output[0].rsplit('.')[0]
    message:
        "[bustools] Generate spliced count matrix from BUS file"
    conda:
        "../envs/bustools.yaml"
    shell:
        "bustools count -o {params.out} -g {input.tsv} -e {input.mat} -t {input.txt} --genecounts {input.bus} 2> {log}"

rule bustools_count_unspliced:
    input:
        tsv = expand("results/busparse/get_velocity_files/{genome}/tr2g.tsv", genome = config["genome"]),
        mat = "results/kallisto/bus/{sample}/matrix.ec",
        txt = "results/kallisto/bus/{sample}/transcripts.txt",
        bus = "results/bustools/capture/{sample}/unspliced.bus"
    output:
        ext = multiext("results/bustools/count/{sample}/unspliced", ".barcodes.txt", ".genes.txt", ".mtx")
    log:
        log = "results/bustools/count/{sample}/unspliced.log"
    params:
        out = lambda wildcards, output: output[0].rsplit('.')[0]
    message:
        "[bustools] Generate unspliced count matrix from BUS file"
    conda:
        "../envs/bustools.yaml"
    shell:
        "bustools count -o {params.out} -g {input.tsv} -e {input.mat} -t {input.txt} --genecounts {input.bus} 2> {log}"
