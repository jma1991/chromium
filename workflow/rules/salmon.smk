# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

salmonAlevin = {
  'A': "--dropseq",
  'B': "--chromiumV3",
  'C': "--chromium",
  'D': "--gemcode",
  'E': "--citeseq",
  'F': "--celseq",
  'G': "--celseq2",
  'H': "--quartzseq2"      
}

rule salmon_index:
    input:
        fas = "results/genomepy/{genome}/{genome}.gentrome.fa",
        txt = "results/genomepy/{genome}/{genome}.decoys.txt"
    output:
        idx = directory("results/salmon/index/{genome}")
    threads:
        16
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon index -t {input.fas} -i {output.idx} -p {threads} -d {input.txt} --keepDuplicates"

rule salmon_alevin:
    input:
        idx = "results/salmon/index/GRCm38.p6",
        fq1 = lambda wildcards: pep.subsample_table.loc[pep.subsample_table['sample_name'] == wildcards.sample, "read1"],
        fq2 = lambda wildcards: pep.subsample_table.loc[pep.subsample_table['sample_name'] == wildcards.sample, "read2"],
        tsv = "results/eisar/GRCm38.p6/GRCm38.p6.tx2gene.tsv",
        txt = ["results/gffread/GRCm38.p6/GRCm38.p6.mrna.txt", "results/gffread/GRCm38.p6/GRCm38.p6.rrna.txt"]
    output:
        mat = "results/salmon/alevin/{sample}/alevin/quants_mat.gz"
    params:
        out = "results/salmon/alevin/{sample}"
    threads:
        16
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon alevin -l ISR -i {input.idx} -1 {input.fq1} -2 {input.fq2} -o {params.out} -p {threads} --tgMap {input.tsv} --chromiumV3 --mrna {input.txt[0]} --rrna {input.txt[1]} --dumpFeatures"
