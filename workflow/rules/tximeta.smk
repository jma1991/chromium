# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule tximeta:
    input:
        fas = "gencode.vM24.annotation.expanded.fa",
        gtf = "gencode.vM24.annotation.expanded.gtf"
    output:
        ""
    script:
        "../scripts/tximeta.R"