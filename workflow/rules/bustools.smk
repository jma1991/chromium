# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule bustools_correct:
    input:
        txt = "whitelist.txt"
        bus = "results/kallisto/{sample}"
    output:
        bus = "results/bustools/{sample}"
    message:
        "[bustools] Error correct BUS file"
    conda:
        "../rules/bustools.yaml"
    shell:
        "bustools correct -o {output.bus} -w {input.txt} {input.bus}

rule bustools_sort:
    input:
        bus = "results/bustools/{sample}.correct.bus"
    output:
        bus = "results/bustools/{sample}.correct.sort.bus"
    message:
        "[bustools] Sort BUS file by barcode and UMI"
    conda:
        "../rules/bustools.yaml"
    shell:
        "bustools sort -t {threads} -o {output.bus}"

rule bustools_capture:
    input:
        mat = "results/kallisto/{sample}/matrix.ec",
        txt = "results/genomepy/transcripts.txt",
        bus = "results/bustools/{sample}.correct.sort.bus"
    output:
        bus = "results/bustools/{sample}.correct.sort.{capture}.bus"
    message:
        "[bustools] Capture records from BUS file"
    shell:
        "bustools capture -o {output.bus} -c {capture} -e {input.mat} -t {input.txt} {input.bus}"

rule bustools_count:
    input:
        bus = ""
    output:
        dir = directory("results/bustools/count/{sample}")
    message:
        "[bustools] Generate count matrix from BUS file"
    shell:
        "bustools count -o {output.dir} -g {input.txt} -e {input.eqc} -t {input.txt} --genecounts {input.bus}"


