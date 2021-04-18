# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule star_index:
    input:
        fas = "results/genomepy/{genome}/{genome}.fa",
        gtf = "results/genomepy/{genome}/{genome}.annotation.gtf"
    output:
        dir = directory("results/star/index/{genome}")
    params:
        arg = star_index_params
    log:
        log = "results/star/index/{genome}/Log.err"
    message:
        "[star] Build STAR index"
    threads:
        16
    conda:
        "../envs/star.yaml"
    shell:
        "STAR"
        " --runMode genomeGenerate"
        " --runThreadN {threads}"
        " --genomeDir {output.dir}"
        " --genomeFastaFiles {input.fas}"
        " --sjdbGTFfile {input.gtf}"
        " --sjdbOverhang {params.arg[sjdbOverhang]}"
        " 2> {log}"

rule star_align:
    input:
        idx = expand("results/star/index/{genome}", genome = config["genome"]),
        gtf = expand("results/genomepy/{genome}/{genome}.annotation.gtf", genome = config["genome"]),
        fq1 = lambda wildcards: units.loc[wildcards.sample, "read1"],
        fq2 = lambda wildcards: units.loc[wildcards.sample, "read2"],
        txt = expand("resources/barcodes/{chemistry}.txt", chemistry = config["chemistry"])
    output:
        dir = directory("results/star/align/{sample}")
    params:
        arg = star_align_params,
        fq1 = lambda wildcards, input: ",".join(input.fq1),
        fq2 = lambda wildcards, input: ",".join(input.fq2),
        out = "results/star/align/{sample}/"
    log:
        log = "results/star/align/{sample}/Log.err"
    message:
        "[star] Run STAR align"
    threads:
        16
    conda:
        "../envs/star.yaml"
    shell:
        "STAR"
        " --runMode alignReads"
        " --runThreadN {threads}"
        " --genomeDir {input.idx}"
        " --sjdbGTFfile {input.gtf}"
        " --readFilesIn {params.fq2} {params.fq1}"
        " --readFilesCommand gunzip -c"
        " --outFileNamePrefix {params.out}"
        " --outSAMtype BAM SortedByCoordinate"
        " --soloType CB_UMI_Simple"
        " --soloCBwhitelist {input.txt}"
        " --soloFeatures Gene Velocyto"
        " --soloCBstart {params.arg[soloCBstart]}"
        " --soloCBlen {params.arg[soloCBlen]}"
        " --soloUMIstart {params.arg[soloUMIstart]}"
        " --soloUMIlen {params.arg[soloUMIlen]}"
        " --soloBarcodeReadLength {params.arg[soloBarcodeReadLength]}"
        " 2> {log}"