# Chromium <img align="right" width="200" src="images/roundel.png">

A Snakemake workflow to process scRNA-seq data from 10x Genomics

## Contents

* [Description](#overview)
* [Installation](#installation)
* [Usage](#usage)
* [Documentation](#documentation)
* [Contributing](#contributing)
* [Authors](#authors)
* [Tests](#tests)
* [Acknowledgements](#acknowledgements)
* [License](#license)

## Description

Chromium is a Snakemake workflow to pre-process 3' single cell gene expression data from the 10x Genomics platform. It is compatible with 10xv2 and 10xv3 chemistry and features three different methods to obtain spliced and unspliced abundance estimates.


## Installation

Chromium requires the following to be installed:

- [conda](https://docs.conda.io/en/latest/index.html)
- [snakedeploy](https://snakedeploy.readthedocs.io/en/latest/)
- [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

Install Chromium with snakedeploy:

   ```console
   $ snakedeploy deploy-workflow https://github.com/jma1991/chromium path/to/project
   ```

## Usage

To run Chromium, follow these instructions:

1. Create workflow configuration

   ```console
   $ vim config/config.yaml
   ```

2. Create samples table

   ```console
   $ vim config/samples.csv
   ```

3. Create units table

   ```console
   $ vim config/units.csv
   ```

4. Test configuration by performing a dry-run

   ```console
   $ snakemake -n
   ```

5. Execute workflow and deploy software dependencies

    ```console
    $ snakemake --cores all --use-conda
    ```


## Roadmap

Chromium 

- Support 10xv1 chemistry
- 



## Documentation

Full documentation is available [here](workflow/documentation.md)


## Features

Chromium has the following:

- 
- compatible with



## Support

If you need any support, please open an [issue](https://github.com/jma1991/scrnaseq/issues) and apply either the *help wanted* or *question* label.

## Feedback

If you have any feedback, please open an [issue](https://github.com/jma1991/scrnaseq/issues) and apply the appropriate label

## Contributing

To contribute to Chromium, clone this repository locally and commit your code on a separate branch. Please generate unit tests for your code and run the linter before opening a pull request:

```console
$ snakemake --generate-unit-tests   # generate unit tests
$ snakemake --lint                  # run the linter
```

You can find more details in the [Contributing](CONTRIBUTING.md) guide. Participation in this open source project is subject to a [Code of Conduct](CODE_OF_CONDUCT.md).

## Authors

Chromium was developed by [James Ashmore](https://www.github.com/jma1991) but has benefited from contributions by the following:

- [Benjamin Southgate](#)
- [Alastair Kilpatrick](#)

If you would like to be added to this list, please open a [pull request](https://github.com/jma1991/scrnaseq/pulls) with your contribution.

## Tests

Test cases are in the `.test` directory. They are automatically executed via  
continuous integration with GitHub Actions.

## Citation

If you use Chromium in your research, please cite using

## Used By

Chromium is used by the following companies and institutes:

- [The Centre for Regenerative Medicine (The University of Edinburgh)](https://www.ed.ac.uk/regenerative-medicine)

If you would like to be added to this list, please open a [pull request](https://github.com/jma1991/scrnaseq/pulls) with your information.

## Related

Here are some related projects:

- [Hoohm/dropSeqPipe](https://github.com/Hoohm/dropSeqPipe)
- [snakemake-workflows/single-cell-rna-seq](https://github.com/snakemake-workflows/single-cell-rna-seq)
- [crazyhottommy/pyflow-cellranger](https://github.com/crazyhottommy/pyflow-cellranger)

## Acknowledgements

The wokflow was motivated by the following projects:

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
