# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule alevinqc:
    input:
        dir = "results/salmon/alevin/{sample}"
    output:
        pdf = "results/alevinqc/{sample}/alevinReport.pdf"
    script:
        "../scripts/alevinqc.R"