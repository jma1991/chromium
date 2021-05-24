# Author: James Ashmore
# Copyright: Copyright 2020, James Ashmore
# Email: jashmore@ed.ac.uk
# License: MIT

def salmon_alevin_params(wildcards):

    d = {
        "10xv1": {
            "protocol": "--gemcode"
        },
        "10xv2": {
            "protocol": "--chromium"
        },
        "10xv3": {
            "protocol": "--chromiumV3"
        }
    }

    k = config["chemistry"]

    v = d.get(k)

    return v


def star_index_params(wildcards):

    d = {
        "10xv1": {
            "sjdbOverhang": 97
        },
        "10xv2": {
            "sjdbOverhang": 97
        },
        "10xv3": {
            "sjdbOverhang": 90
        }
    }

    k = config["chemistry"]

    v = d.get(k)

    return v


def star_align_params(wildcards):

    d = {
        "10xv1": {
            "sjdbOverhang": 97,
            "soloCBstart": 1,
            "soloCBlen": 14,
            "soloUMIstart": 15,
            "soloUMIlen": 10,
            "soloBarcodeReadLength": 24
        },
        "10xv2": {
            "sjdbOverhang": 97,
            "soloCBstart": 1,
            "soloCBlen": 16,
            "soloUMIstart": 17,
            "soloUMIlen": 10,
            "soloBarcodeReadLength": 26
        },
        "10xv3": {
            "sjdbOverhang": 90,
            "soloCBstart": 1,
            "soloCBlen": 16,
            "soloUMIstart": 17,
            "soloUMIlen": 12,
            "soloBarcodeReadLength": 28
        }
    }

    k = config["chemistry"]

    v = d.get(k)

    return v
