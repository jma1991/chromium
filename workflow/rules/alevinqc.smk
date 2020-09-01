# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule alevinqc:
    input:
        mat = "results/salmon/alevin/{sample}/alevin/quants_mat.gz"
    output:
        html = "results/alevinqc/{sample}.html"
    params:
        dir = "results/salmon/alevin/{sample}"
    message:
        "[alevinQC] Generate alevin summary report"
    conda:
        "../envs/bioconductor-alevinqc.yaml"
    script:
        "../scripts/alevinqc.R"