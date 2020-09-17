# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

rule curl_barcodes:
    output:
        txt = "resources/barcodes/3M-february-2018.txt"
    log:
        out = "resources/barcodes/3M-february-2018.log"
    params:
        url = "https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"
    shell:
        "curl {params.url} | gunzip 1> {output.txt} 2> {log}"
