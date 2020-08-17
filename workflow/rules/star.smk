# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

soloType = {
    'A' = "CBUMISimple",
    'B' = "CBUMIComplex",
    'C' = "CBsamTagOut",
    'D' = "SmartSeq"
}

rule star_index:
    input:
        ""
    output:
        ""
    shell:
        ""

rule star_solo:
    input:
        ""
    output:
        "features.tsv", "barcodes.tsv", "matrix.mtx"
    shell:
        "STAR --runMode alignReads --runThreadN {threads} --genomeDir {input.idx} --sjdbGTFfile {input.gff} --readFilesIn {input.fq1} {input.fq2} --readFilesCommand zcat --outFileNamePrefix {params.out} --outSAMtype BAM SortedByCoordinate --twopassMode Basic --soloType Droplet --soloCBwhitelist {whitelist.txt}"
