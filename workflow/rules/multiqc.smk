# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule multiqc:
    input:
        "results/kallisto", "results/salmon", "results/star"
    output:
        "results/multiqc"
    conda:
        "../rules/multiqc.yaml"
    shell:
        "multiqc"