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

rule salmon_gentrome:
    input:
        "results/eisar/{genome}/{genome}.annotation.expanded.fa",
        "results/genomepy/{genome}/{genome}.fa"
    output:
        "results/salmon/genome/{genome}/gentrome.fa"
    shell:
        "cat {input} > {output}"

rule salmon_decoys:
	input:
		"results/genomepy/{genome}/{genome}.fa.fai"
	output:
		"results/salmon/genome/{genome}/decoys.txt"
	shell:
		"cut -f 1 {input} > {output}"

rule salmon_mrna:
	input:
		"results/genomepy/{genome}/{genome}.annotation.gtf"
	output:
		"results/salmon/genome/{genome}/mrna.txt"
	shell:
		"gffread {input} --table @chr,gene_id | awk '$1 == MT' | cut -f 2 | sort -u > {output}"

rule salmon_rrna:
	input:
		"results/genomepy/{genome}/{genome}.annotation.gtf"
	output:
		"results/salmon/genome/{genome}/rrna.txt"
	shell:
		"gffread {input} --table gene_biotype,gene_id | awk '$1 == rRNA' | cut -f 2 | sort -u > {output}"

rule salmon_index:
    input:
        fas = "results/salmon/genome/{genome}/gentrome.fa",
        txt = "results/salmon/genome/{genome}/decoys.txt"
    output:
        idx = directory("results/salmon/index/{genome}")
    shell:
        "salmon index -t {input.fas} -i {output.idx} -p {threads} -d {input.txt}"

rule salmon_alevin:
    input:
        idx = "results/salmon/index/GRCm38.p6",
        fq1 = lambda wildcards: units.loc[units['sample'] == wildcards.sample, "fq1"],
        fq2 = lambda wildcards: units.loc[units['sample'] == wildcards.sample, "fq2"],
        tsv = "results/eisar/GRCm38.p6/GRCm38.p6.annotation.expanded.tx2gene.tsv",
        txt = ["results/salmon/genome/GRCm38.p6/mrna.txt", "results/salmon/genome/GRCm38.p6/rrna.txt"]
    output:
        dir = directory("results/salmon/alevin/{sample}")
    threads:
        16
    shell:
        "salmon alevin -l ISR -i {input.idx} -1 {input.fq1} -2 {input.fq2} -o {output.dir} -p {threads} --tgMap {input.tsv} --chromiumV3 --mrna {input.txt[0]} --rrna {input.txt[1]}"
