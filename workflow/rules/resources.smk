# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule barcodes_10xv1:
    output:
        txt = "resources/barcodes/10xv1.txt"
    params:
        url = "https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/737K-april-2014_rc.txt"
    log:
        err = "resources/barcodes/10xv1.log"
    message:
        "[resources] Download 10xv1 barcodes"
    conda:
        "../envs/curl.yaml"
    shell:
        "curl {params.url} > {output.txt} 2> {log.err}"

rule barcodes_10xv2:
    output:
        txt = "resources/barcodes/10xv2.txt"
    params:
        url = "https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/737K-august-2016.txt"
    log:
        err = "resources/barcodes/10xv2.log"
    message:
        "[resources] Download 10xv2 barcodes"
    conda:
        "../envs/curl.yaml"
    shell:
        "curl {params.url} > {output.txt} 2> {log.err}"

rule barcodes_10xv3:
    output:
        txt = "resources/barcodes/10xv3.txt"
    params:
        url = "https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"
    log:
        err = "resources/barcodes/10xv3.log"
    message:
        "[resources] Download 10xv3 barcodes"
    conda:
        "../envs/curl.yaml"
    shell:
        "curl {params.url} | gunzip -c > {output.txt} 2> {log.err}"
