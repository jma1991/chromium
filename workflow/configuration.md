## Configuration <img align="right" width="200" src="../images/roundel.png">

### Contents

* [Overview](#overview)
* [Create a workflow config](#create-a-config-file)
* [Create a samples table](#create-a-samples-table)
* [Create a units table](#create-a-units-table)

### Overview

Chromium is configured by editing the files in the `config` directory:

- the `config.yaml` file contains the workflow settings
- the `samples.csv` file contains the samples table
- the `units.csv` file contains the units table

The workflow will throw an error if these files are missing or do not contain  
the required information.

### Create a workflow config

The workflow config is a YAML file containing information about the  
workflow parameters:

* each line contains a name:value pair
* each name:value pair corresponds to a workflow parameter

The workflow config must contain the following:

| Name | Description | Example |
| --- | --- | --- |
| samples | Path to samples.csv | config/samples.csv |
| units | Path to units.csv | config/units.csv |
| source | Transcriptome source | Ensembl |
| organism | Species name | Homo sapiens |
| release | Release number | 101 |
| genome | Genome assembly | GRCh38.p13 |
| chemistry | Chemistry version | 10xv3 |

Example of a valid workflow config:

```yaml
samples: "config/samples.csv"
units: "config/units.csv"
source: "Ensembl"
organism: "Homo sapiens"
release: "101"
genome: "GRCh38.p13"
chemistry: "10xv3"
```

### Create a samples table

The samples table is a CSV file containing information about the samples in  
your data set:

* each row corresponds to one sample
* each column corresponds to one attribute

For each sample, you must provide the following:

| Column | Description | Example |
| --- | --- | --- |
| sample | Sample name | S1 |

Example of a valid samples table:

```
sample,condition
S1,control
S2,treatment
```

### Create a units table

The units table is a CSV file containing information about the sequencing units  
in your data set:

* each row corresponds to one unit
* each column corresponds to one attribute

For each unit, you must provide the following:

| Column | Description | Example |
| --- | --- | ---|
| sample | Sample name | S1 |
| unit | Unit name | L001 |
| read1 | Read 1 file | S1_L001_1.fastq.gz |
| read2 | Read 2 file | S1_L001_2.fastq.gz |

Example of valid units table:

```
sample,unit,read1,read2
S1,L001,S1_L001_1.fastq.gz,S1_L001_2.fastq.gz
S1,L002,S1_L002_1.fastq.gz,S1_L002_2.fastq.gz
S2,L001,S2_L001_1.fastq.gz,S2_L001_2.fastq.gz
S2,L002,S2_L002_1.fastq.gz,S2_L002_2.fastq.gz
```