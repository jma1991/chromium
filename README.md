# Chromium <img align="right" width="200" src="images/roundel.png">

Processing of scRNA-seq data

[![Snakemake](https://img.shields.io/badge/snakemake-6.2.1-brightgreen.svg)](https://snakemake.readthedocs.io)
[![MIT](https://img.shields.io/badge/license-MIT-blue)](https://opensource.org/licenses/MIT)

## Contents

* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
* [Acknowledgements](#acknowledgements)
* [Authors](#authors)
* [Contributing](#contributing)
* [Documentation](#documentation)
* [FAQ](#faq)
* [Features](#features)
* [Installation](#installation)
* [License](#license)
* [Tests](#tests)
* [Citation](#citation)

## Introduction

Chromium is a Snakemake workflow to pre-process 3' single cell gene expression data from the 10x Genomics platform. It features three different quantification methods to obtain both spliced and unspliced abundance estimates summarised by gene:

1. kallisto|bustools
2. Alevin
3. STARsolo

## Installation

Chromium requires the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system and [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) package management system to be installed. For instructions, please read their respective documentation. Once available, Chromium and all of its dependencies can be installed with the following commands:

1. Clone workflow into working directory: 
    
    ```console
    $ git clone https://github.com/jashmore/chromium.git
    ```

2. Change to workflow directory:

    ```console
    $ cd chromium
    ```

3. Install required conda environments:

    ```console
    $ snakemake --use-conda --conda-create-envs-only
    ```

## Usage

1. Clone workflow into working directory:

    ```console
    $ git clone https://github.com/jashmore/chromium.git
    ```

2. Change to workflow directory:

    ```console
    $ cd chromium
    ```

3. Edit workflow configuration:

    ```console
    $ vim config/config.yaml   # use your favourite text editor
    ```

4. Edit samples table:

    ```console
    $ vim config/samples.csv   # use your favourite text editor
    ```

5. Edit units table:

    ```console
    $ vim config/units.csv   # use your favourite text editor
    ```

6. Test configuration by performing a dry-run:

    ```console
    $ snakemake -n
    ```

7. Execute workflow and deploy software dependencies via conda:

    ```console
    $ snakemake --use-conda
    ```

## Acknowledgements

- [RNA velocity with kallisto | bus and velocyto.R](https://bustools.github.io/BUS_notebooks_R/velocity.html)
- [Alevin velocity](https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/)
- [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [Writing a Friendly README](https://rowanmanning.com/posts/writing-a-friendly-readme/)
- [readme.so](https://readme.so)

## Documentation

See [`workflow/documentation.md`](workflow/documentation.md) for the full project documentation.

## FAQ

#### Which quantification workflow should I use?

I would suggest there is no best workflow, each one captures a unique aspect of the data by the counting strategies they have implemented. For a more in-depth discussion, please refer to this research article by Soneson and colleagues: https://doi.org/10.1371/journal.pcbi.1008585

#### What reference genomes are supported?

The workflow uses `genomepy` and `gffread` to download and parse the user-specified reference genome and annotation. Theoretically, any genome release compatible with these software should be supported.

#### How do I combine multiple sequencing runs?

If the sequencing runs were performed across multiple lanes on the same date, it is unlikely that a batch effect is present and I would recommend quantifying the files all together. Below is an example of how to specify multiple sequencing runs jointly for a given sample: 

```
sample,unit,read1,read2
S1,L001,S1_L001.fastq.gz,S1_L001.fastq.gz
S1,L002,S1_L002.fastq.gz,S1_L002.fastq.gz
```

Alternatively, if the sequencing runs were performed on different machines and different dates, there is potential for a batch effect and I would recommend quantifying the files separately until this can be investigated. Below is an example of how to specify multiple sequencing runs independently for a given sample:

```
sample,unit,read1,read2
S1_L001,L001,S1_L001.fastq.gz,S1_L001.fastq.gz
S1_L002,L002,S1_L002.fastq.gz,S1_L002.fastq.gz
```

#### 






## Authors

Chromium was initially developed by [James Ashmore](https://www.github.com/jma1991) but has benefited from contributions by many individuals in the community:

- [Benjamin Southgate](#)
- [Alastair Kilpatrick](#)


## Features

Chromium offers the following features:

- Performs spliced and unspliced transcript quantification
- Implements multiple pre-processing workflows: Kallisto-bustools, Salmon-alevin, and STAR-solo
- Supports Single Cell 3' v1, v2, v3, and LT chemistry
- Outputs SingleCellExperiment objects for downstream analysis

## Contributing

To contribute to Chromium, clone this repository locally and commit your code on a separate branch. Please generate unit tests for your code and run the linter before opening a pull-request:

```console
$ snakemake --generate-unit-tests   # generate unit tests
$ snakemake --lint                  # run the linter
```

You can find more details in our [Contributing](CONTRIBUTING.md) guide. Participation in this open source project is subject to a [Code of Conduct](CODE_OF_CONDUCT.md).

## Tests

Test cases are in the `.test` directory. They are automatically executed via  
continuous integration with GitHub Actions.

## Citation

If you use Chromium in your research, please cite using

## License

Chromium is licensed under the [MIT](LICENSE) license.  
Copyright &copy; 2020, James Ashmore