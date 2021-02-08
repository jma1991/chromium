# SCRNA-SEQ <img align="right" width="200" src="images/roundel.png">

A workflow for scRNA-seq analysis in Snakemake.

[![Snakemake][shield-snakemake]](https://snakemake.readthedocs.io)
[![MIT license][shield-license]](https://choosealicense.com/licenses/mit)

Table of Contents
-----------------

  * [Introduction](#introduction)
  * [Requirements](#requirements)
  * [Usage](#usage)
  * [Contributing](#contributing)
  * [Thanks](#thanks)
  * [License](#license)

Introduction
------------

This workflow is a bioinformatics analysis pipeline for single-cell RNA sequencing data. The workflow is built using [Snakemake - a scalabale bioinformatics workflow engine](https://doi.org/10.1093/bioinformatics/bts480)


Requirements
------------

This workflow requires the following software to run:

  * [Snakemake][snakemake]
  * [Conda][code]

Usage
-----


Ripple is easiest to use when installed with [npm][npm]:


Clone workflow into working directory:

```sh
git clone https://github.com/jma1991/scrnaseq.git
```

Execute workflow and deploy software dependencies via conda:

```sh
snakemake --use-conda
```

Configuration
-------------

Configure the workflow by editing the files in the `config` directory:

- `config.yaml` is a YAML file containing the workflow metadata.

- `samples.csv` is a CSV file containing the sample metadata.

- `units.csv` is a CSV file contains the unit metadata.

Contributing
------------

To contribute to the workflow, clone this repository locally and commit your code on a separate branch. Please generate unit tests for your code, and run the linter before opening a pull-request:

```sh
snakemake --generate-unit-tests # generate unit tests
snakemake --lint # run the linter
```

You can find more detail in our [Contributing Guide](#). Participation in this open source project is subject to a [Code of Conduct](#).

Thanks
------

I would like to thank Johannes Köster for developing the Snakemake workflow engine and Istvan Albert for writing the biostar handbook.

License
-------

This workflow is licensed under the [MIT](#) license.  
Copyright &copy; 2020, James Ashmore


[shield-snakemake]: https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg
[shield-license]: https://img.shields.io/badge/license-MIT-blue.svg
