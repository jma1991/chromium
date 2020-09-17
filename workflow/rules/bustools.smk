# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule bustools_correct:
    input:
        txt = "resources/barcodes/3M-february-2018.txt",
        bus = "results/kallisto/{sample}/output.bus"
    output:
        bus = "results/bustools/{sample}/output.correct.bus"
    log:
        log = "results/bustools/{sample}/output.correct.log"
    message:
        "[bustools] Error correct BUS file"
    conda:
        "../envs/bustools.yaml"
    shell:
        "bustools correct -o {output.bus} -w {input.txt} {input.bus} 2> {log}"

rule bustools_sort:
    input:
        bus = "results/bustools/{sample}/output.correct.bus"
    output:
        bus = "results/bustools/{sample}/output.correct.sort.bus"
    log:
        log = "results/bustools/{sample}/output.correct.sort.log"
    message:
        "[bustools] Sort BUS file by barcode and UMI"
    conda:
        "../envs/bustools.yaml"
    shell:
        "bustools sort -t {threads} -o {output.bus} {input.bus} 2> {log}"

rule bustools_capture_spliced:
    input:
        txt = "results/busparse/GRCm38.p6/introns_tx_to_capture.txt",
        mat = "results/kallisto/{sample}/matrix.ec",
        txi = "results/kallisto/{sample}/transcripts.txt",
        bus = "results/bustools/{sample}/output.correct.sort.bus"
    output:
        bus = "results/bustools/{sample}/output.correct.sort.spliced.bus"
    log:
        log = "results/bustools/{sample}/output.correct.sort.spliced.log"
    message:
        "[bustools] Capture spliced reads from BUS file"
    conda:
        "../envs/bustools.yaml"
    shell:
        "bustools capture -s -x -o {output.bus} -c {input.txt} -e {input.mat} -t {input.txi} {input.bus} 2> {log}"

rule bustools_capture_unspliced:
    input:
        txt = "results/busparse/GRCm38.p6/cDNA_tx_to_capture.txt",
        mat = "results/kallisto/{sample}/matrix.ec",
        txi = "results/kallisto/{sample}/transcripts.txt",
        bus = "results/bustools/{sample}/output.correct.sort.bus"
    output:
        bus = "results/bustools/{sample}/output.correct.sort.unspliced.bus"
    log:
        log = "results/bustools/{sample}/output.correct.sort.unspliced.log"
    message:
        "[bustools] Capture spliced reads from BUS file"
    conda:
        "../envs/bustools.yaml"
    shell:
        "bustools capture -s -x -o {output.bus} -c {input.txt} -e {input.mat} -t {input.txi} {input.bus} 2> {log}"

rule bustools_count_spliced:
    input:
        tsv = "results/busparse/GRCm38.p6/tr2g.tsv",
        mat = "results/kallisto/{sample}/matrix.ec",
        txt = "results/kallisto/{sample}/transcripts.txt",
        bus = "results/bustools/{sample}/output.correct.sort.spliced.bus"
    output:
        ext = multiext("results/bustools/{sample}/output.correct.sort.spliced", ".barcodes.txt", ".genes.txt", ".mtx")
    log:
        log = "results/bustools/{sample}/output.correct.sort.spliced.log"
    params:
        out = "results/bustools/{sample}/output.correct.sort.spliced"
    message:
        "[bustools] Generate spliced count matrix from BUS file"
    conda:
        "../envs/bustools.yaml"
    shell:
        "bustools count -o {params.out} -g {input.tsv} -e {input.mat} -t {input.txt} --genecounts {input.bus} 2> {log}"

rule bustools_count_unspliced:
    input:
        tsv = "results/busparse/GRCm38.p6/tr2g.tsv",
        mat = "results/kallisto/{sample}/matrix.ec",
        txt = "results/kallisto/{sample}/transcripts.txt",
        bus = "results/bustools/{sample}/output.correct.sort.unspliced.bus"
    output:
        ext = multiext("results/bustools/{sample}/output.correct.sort.unspliced", ".barcodes.txt", ".genes.txt", ".mtx")
    log:
        log = "results/bustools/{sample}/output.correct.sort.unspliced.log"
    params:
        out = "results/bustools/{sample}/output.correct.sort.unspliced"
    message:
        "[bustools] Generate unspliced count matrix from BUS file"
    conda:
        "../envs/bustools.yaml"
    shell:
        "bustools count -o {params.out} -g {input.tsv} -e {input.mat} -t {input.txt} --genecounts {input.bus} 2> {log}"
