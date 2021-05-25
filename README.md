# Chromium <img align="right" width="200" src="images/roundel.png">

A Snakemake workflow to process scRNA-seq data from 10x Genomics

[![DOI](https://zenodo.org/badge/287750880.svg)](https://zenodo.org/badge/latestdoi/287750880)

## Contents

* [Overview](#overview)
* [Installation](#installation)
* [Usage](#usage)
* [Documentation](#documentation)
* [Contributing](#contributing)
* [Authors](#authors)
* [Tests](#tests)
* [Acknowledgements](#acknowledgements)
* [License](#license)

## Overview

Chromium is a Snakemake workflow to process 3' single cell RNA sequencing data from the 10x Genomics platform. It is compatible with 10xv2 and 10xv3 chemistry and features three different quantification methods to obtain both spliced and unspliced abundance estimates:

* [Kallisto/Bustools](https://doi.org/10.1038/s41587-021-00870-2)
* [Alevin](https://doi.org/10.1186/s13059-019-1670-y)
* [STARsolo](https://doi.org/10.1101/2021.05.05.442755)

## Installation

Chromium and all of its dependencies can be installed via the [mamba](https://github.com/mamba-org/mamba) package manager:

1. Install Snakemake and Snakedeploy

   ```console
   $ mamba create -c bioconda -c conda-forge --name snakemake snakemake snakedeploy
   ```

2. Activate the Snakemake environment

   ```console
   $ mamba activate snakemake
   ```

3. Create a project directory

   ```console
   $ mkdir -p path/to/project
   ```

4. Deploy the workflow in the project directory

   ```console
   $ snakedeploy deploy-workflow https://github.com/snakemake-workflows/chromium path/to/project
   ```

## Usage

1. Create workflow configuration

   ```console
   $ vim config/config.yaml   # workflow parameters
   ```

2. Create samples table

   ```console
   $ vim config/samples.csv   # sample metadata
   ```

3. Create units table

   ```console
   $ vim config/units.csv   # unit metdata
   ```

4. Test configuration by performing a dry-run

   ```console
   $ snakemake -n
   ```

5. Execute workflow and deploy software dependencies

    ```console
    $ snakemake --cores all --use-conda
    ```

*For more information, see the [Usage](workflow/documentation.md#usage) section of the documentation.*

## Documentation

Full documentation for Chromium is available [here](workflow/documentation.md)

## Support

If you need help, open an [issue](https://github.com/jma1991/scrnaseq/issues) with one of the following labels:

- help wanted (extra attention is needed)
- question (further information is requested)

## Feedback

If you have any suggestions, open an [issue](https://github.com/jma1991/scrnaseq/issues) with one of the following labels:

- documentation (improvements or additions to documentation)
- enhancement (new feature or request)

## Contributing

To contribute to Chromium, clone this repository locally and commit your code on a separate branch. Please generate unit tests for your code and run the linter before opening a pull request:

```console
$ snakemake --generate-unit-tests   # generate unit tests
$ snakemake --lint                  # run the linter
```

You can find more details in the [Contributing](CONTRIBUTING.md) guide. 

Participation in this project is subject to a [Code of Conduct](CODE_OF_CONDUCT.md).

## Authors

Chromium was developed by [James Ashmore](https://www.github.com/jma1991) but has benefited from contributions by the following:

- [Benjamin Southgate](#)
- [Alastair Kilpatrick](#)

If you would like to be added to this list, please open a [pull request](https://github.com/jma1991/scrnaseq/pulls) with your contribution.

## Citation

If you use Chromium in your research, please cite using the DOI: [10.5281/zenodo.4783308](https://doi.org/10.5281/zenodo.4783309)

## Used By

Chromium is used by the following companies and institutes:

- [The Centre for Regenerative Medicine (The University of Edinburgh)](https://www.ed.ac.uk/regenerative-medicine)

If you would like to be added to this list, please open a [pull request](https://github.com/jma1991/scrnaseq/pulls) with your information.

## Acknowledgements

The methods were chosen based on this research paper:

> Soneson C, Srivastava A, Patro R, Stadler MB (2021) Preprocessing choices affect RNA velocity results for droplet scRNA-seq data. PLoS Comput Biol 17(1): e1008585. https://doi.org/10.1371/journal.pcbi.1008585

The workflow was motivated by the following projects:

- [nf-core/scrnaseq](https://github.com/nf-core/scrnaseq)
- [maxplanck-ie/snakepipes](https://github.com/maxplanck-ie/snakepipes)
- [10XGenomics/cellranger](https://github.com/10XGenomics/cellranger)

The documentation was informed by the following articles:

- [easiest way to create a readme](https://readme.so)
- [writing a friendly readme](https://rowanmanning.com/posts/writing-a-friendly-readme/)
- [writing well for the web](https://www.gov.uk/guidance/content-design/writing-for-gov-uk)

## License

Chromium is licensed under the [MIT](LICENSE.md) license.  
Copyright &copy; 2020, James Ashmore
