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
    threads:
        16
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.dir} --genomeFastaFiles {input.fas} --sjdbGTFfile {input.gtf} --sjdbOverhang 100"

rule star_solo:
    input:
        idx = "results/star/index/GRCm38.p6",
        gtf = "results/genomepy/GRCm38.p6/GRCm38.p6.annotation.gtf",
        fq1 = lambda wildcards: units.loc[units["sample"] == wildcards.sample, "fq1"],
        fq2 = lambda wildcards: units.loc[units["sample"] == wildcards.sample, "fq2"],
        txt = "workflow/resources/barcodes/3M-february-2018.txt"
    output:
        dir = directory("results/star/solo/{sample}")
    params:
        fq1 = lambda wildcards, input: ",".join(input.fq1),
        fq2 = lambda wildcards, input: ",".join(input.fq2),
        out = "results/star/solo/{sample}/"
    threads:
        16
    shell:
        "STAR --runMode alignReads --runThreadN {threads} --genomeDir {input.idx} --sjdbGTFfile {input.gtf} --readFilesIn {params.fq2} {params.fq1} --readFilesCommand gunzip -c --outFileNamePrefix {params.out} --outSAMtype BAM SortedByCoordinate --twopassMode Basic --soloType CB_UMI_Simple --soloCBwhitelist {input.txt} --soloFeatures Gene Velocyto --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 151"
