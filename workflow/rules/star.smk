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
    log:
        log = "results/star/index/{genome}/.snakemake.log"
    threads:
        16
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.dir} --genomeFastaFiles {input.fas} --sjdbGTFfile {input.gtf} --sjdbOverhang 150 2> {log}"

rule star_solo:
    input:
        idx = "results/star/index/GRCm38.p6",
        gtf = "results/genomepy/GRCm38.p6/GRCm38.p6.annotation.gtf",
        fq1 = lambda wildcards: pep.subsample_table.loc[pep.subsample_table['sample_name'] == wildcards.sample, "read1"],
        fq2 = lambda wildcards: pep.subsample_table.loc[pep.subsample_table['sample_name'] == wildcards.sample, "read2"],
        txt = "resources/barcodes/3M-february-2018.txt"
    output:
        bam = "results/star/solo/{sample}/Aligned.sortedByCoord.out.bam"
    log:
        log = "results/star/solo/{sample}/.snakemake.log"
    params:
        fq1 = lambda wildcards, input: ",".join(input.fq1),
        fq2 = lambda wildcards, input: ",".join(input.fq2),
        out = "results/star/solo/{sample}/"
    threads:
        16
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode alignReads --runThreadN {threads} --genomeDir {input.idx} --sjdbGTFfile {input.gtf} --readFilesIn {params.fq2} {params.fq1} --readFilesCommand gunzip -c --outFileNamePrefix {params.out} --outSAMtype BAM SortedByCoordinate --soloType CB_UMI_Simple --soloCBwhitelist {input.txt} --soloFeatures Gene Velocyto --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 151 2> {log}"
