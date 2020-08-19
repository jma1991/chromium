# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule txome:
    input:
        fas = "gencode.vM24.annotation.expanded.fa",
        gtf = "gencode.vM24.annotation.expanded.gtf"
    output:
        ""
    script:
        "../scripts/txome.R"



tximeta::makeLinkedTxome(
  indexDir = "gencode.vM24.annotation.expanded.sidx", 
  source = "GENCODE", genome = "GRCm38", 
  organism = "Mus musculus", release = "M24", 
  fasta = "gencode.vM24.annotation.expanded.fa", 
  gtf = "gencode.vM24.annotation.expanded.gtf", 
  write = TRUE, jsonFile = "gencode.vM24.annotation.expanded.json"
)