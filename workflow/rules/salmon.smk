# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule salmon_index:
    input:
        fas = "results/genomepy/{genome}/{genome}.gentrome.fa",
        txt = "results/genomepy/{genome}/{genome}.decoys.txt"
    output:
        idx = directory("results/salmon/index/{genome}")
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
        "salmon index -t {input.fas} -i {output.idx} -p {threads} -d {input.txt} --keepDuplicates 1> {log.out} 2> {log.err}"

rule salmon_alevin:
    input:
        idx = "results/salmon/index/GRCm38.p6",
        fq1 = lambda wildcards: pep.subsample_table.loc[pep.subsample_table['sample_name'] == wildcards.sample, "read1"],
        fq2 = lambda wildcards: pep.subsample_table.loc[pep.subsample_table['sample_name'] == wildcards.sample, "read2"],
        tsv = "results/eisar/GRCm38.p6/GRCm38.p6.tx2gene.tsv",
        txt = ["results/gffread/GRCm38.p6/GRCm38.p6.mrna.txt", "results/gffread/GRCm38.p6/GRCm38.p6.rrna.txt"]
    output:
        mat = "results/salmon/alevin/{sample}/alevin/quants_mat.gz"
    log:
        out = "results/salmon/alevin/{sample}/log.out",
        err = "results/salmon/alevin/{sample}/log.err"
    params:
        out = "results/salmon/alevin/{sample}"
    message:
        "[salmon]"
    threads:
        16
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon alevin -l ISR -i {input.idx} -1 {input.fq1} -2 {input.fq2} -o {params.out} -p {threads} --tgMap {input.tsv} --chromiumV3 --mrna {input.txt[0]} --rrna {input.txt[1]} --dumpFeatures 1> {log.out} 2> {log.err}"