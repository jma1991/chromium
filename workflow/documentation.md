# Documentation <img align="right" width="160" src="../images/roundel.png">

Welcome to the Chromium documentation

## Contents

* [Introduction](#overview)
* [Usage](#execution)
* [Configuration](#configuration)
* [Output](#results)
* [Tests](#test)
* [FAQ](#faq)
* [References](#references)

## Introduction

Chromium is a Snakemake workflow to process single cell gene expression data
from the 10x Genomics platfom.

## Usage

The workflow can be executed using the following command:

```console
$ snakemake --cores all --use-conda
```

This will use all available cores and deploy software dependencies via the conda
package manager. For further information, please refer to the official
[Snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
documentation.

### Update the workflow

Use git to pull the latest release:

   ```console
   $ git pull
   ```

## Configuration

Chromium is configured by editing the files in the config directory:

- config.yaml 
- samples.csv
- units.csv

An error will be thrown if these files are missing or do not contain the
required information.

### Workflow config

The workflow config is a YAML file containing information about the workflow
parameters:

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

### Samples table

The samples table is a CSV file containing information about the biological
samples in your experiment:

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

### Units table

The units table is a CSV file containing information about the sequencing units
in your experiment:

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

## Output

Chromium processes the data using three different quantification methods:

* kallisto|bustools
* Alevin
* STARsolo

For each method, all output files are saved to disk and a SingleCellExperiment  
object containing spliced and unspliced count matrices is generated.

### Output directory

Chromium saves all output files in the `output` directory. The files are organised  
by software and labelled by sample name and data type:

```console
output
│
├── busparse
│   └── {genome}.cDNA_introns.fa
│   ├── {genome}.cDNA_tx_to_capture.txt
│   ├── {genome}.introns_tx_to_capture.txt
│   └── {genome}.tr2g.tsv
│
├── bustools
│   ├── {sample}.correct.bus
│   ├── {sample}.sort.bus
│   ├── {sample}.spliced.bus
│   ├── {sample}.unspliced.bus
│   ├── {sample}.spliced.mtx
│   ├── {sample}.spliced.barcodes.txt
│   ├── {sample}.spliced.genes.txt
│   ├── {sample}.unspliced.mtx
│   ├── {sample}.unspliced.barcodes.txt
│   └── {sample}.unspliced.genes.txt
│
├── eisar
│   ├── {genome}.annotation.gtf
│   ├── {genome}.fa
│   ├── {genome}.features.tsv
│   ├── {genome}.ranges.rds
│   └── {genome}.tx2gene.tsv
│
├── genomepy
│   └── index
│       └── {genome}.idx

│   ├── {genome}.annotation.bed.gz
│   ├── {genome}.annotation.gtf.gz
│   ├── {genome}.fa
│   ├── {genome}.fa.fai
│   ├── {genome}.fa.sizes
│   └── {genome}.gaps.bed
│
├── gffread
│   ├── {genome}.id2name.tsv
│   ├── {genome}.mrna.txt
│   ├── {genome}.rrna.txt
│   └── {genome}.tx2gene.tsv
│
├── kallisto
│   ├── bus
│   │   └── {sample}
│   └── index
│       └── {genome}.idx
│
├── salmon
│   ├── alevin
│   │   └── {sample}
│   ├── genome
│   │   └── {genome}
│   └── index
│       └── {genome}
│
├── singlecellexperiment
│   ├── {sample}.kallisto.rds
│   ├── {sample}.salmon.rds
│   └── {sample}.star.rds
│
└── star
    ├── align
    │   └── {sample}
    └── index
        └── {genome}
```

### Output files

For each rule


#### BUSpaRse

The `busparse` directory contains output files generated by the
`get_velocity_files` and `read_velocity_files` functions from the
[BUSpaRse](https://bioconductor.org/packages/BUSpaRse/) software:

| File | Format | Description |
| --- | --- | --- |
| `{genome}.cDNA_introns.fa` | FASTA | Spliced transcript and intron sequences |
| `{genome}.cDNA_tx_to_capture.txt` | TXT | Transcript IDs of spliced transcripts |
| `{genome}.introns_tx_to_capture.txt` | TXT | Transcript IDs of introns |
| `{genome}.tr2g.tsv` | TSV | Transcript to gene ID table |

#### BUStools

The `bustools` directory contains output files generated by the `correct`, `sort`, `capture`, and `count` commands of the <ins>[BUStools](https://doi.org/10.1093/bioinformatics/btz279)</ins> software:

| File | Format | Description |
| --- | --- | --- |
| `{sample}.correct.bus` | BUS | Corrected BUS file |
| `{sample}.sort.bus` | BUS | Sorted BUS file |
| `{sample}.spliced.bus` | BUS | Spliced BUS file |
| `{sample}.unspliced.bus` | BUS | Unspliced BUS file |
| `{sample}.spliced.barcodes.txt` | TXT | Spliced barcodes |
| `{sample}.spliced.genes.txt` | TXT | Spliced genes |
| `{sample}.spliced.mtx` | MTX | Spliced counts |
| `{sample}.unspliced.barcodes.txt` | TXT | Unspliced barcodes |
| `{sample}.unspliced.genes.txt` | TXT | Unspliced genes |
| `{sample}.unspliced.mtx` | MTX | Unspliced counts |

#### eisaR

The `eisar` directory contains output files generated by the `exportToGtf`,  
`getFeatureRanges`, and `getTx2Gene` functions from the [eisaR](https://bioconductor.org/packages/eisaR/) software:

| File | Format | Description |
| --- | --- | --- |
| `{genome}.annotation.gtf` | GTF | Expanded gene annotation |
| `{genome}.fa` | FASTA | Spliced transcript and intron sequences |
| `{genome}.features.tsv` | TSV | Spliced transcript and intron names |
| `{genome}.ranges.rds` | RDS | Spliced transcript and intron ranges |
| `{genome}.tx2gene.tsv` | TSV | Transcript to gene table |

#### genomepy

The `genomepy` directory contains output files generated by the `install`  
command from the [genomepy](https://github.com/vanheeringen-lab/genomepy) software.

| File | Format | Description |
| --- | --- | --- |
| `{genome}.annotation.bed.gz` | BED | Gene annotation |
| `{genome}.annotation.gtf.gz` | GTF | Gene annotation |
| `{genome}.fa` | FASTA | Genome sequence |
| `{genome}.fa.sizes` | TSV | Chromosome size |
| `{genome}.gaps.bed` | BED | Gap location |
| `README.txt` | TXT | README |

#### GffRead

The `gffread` directory contains output files generated by the <ins>[GffRead](https://f1000research.com/articles/9-304/v2)</ins> software:

| File | Format | Description |
| --- | --- | --- |
| `{genome}.id2name.tsv` | TSV | The gene_id to gene_name annotation table |
| `{genome}.mrna.txt` | TXT | The mRNA gene_id annotation table |
| `{genome}.rrna.txt` | TXT | The rRNA gene_id annotation table |
| `{genome}.tx2gene.tsv` | TSV | The transcript_id to gene_id annotation table |

#### Kallisto

The `kallisto` directory contains output files generated by the `index` and `bus` commands of the <ins>[Kallisto]()</ins> software:

| File | Format | Description |
| --- | --- | --- |
| `bus/{sample}/matrix.ec` | MTX | Equivalence class |
| `bus/{sample}/output.bus` | BUS | Output |
| `bus/{sample}/run_info.json` | JSON | Run information |
| `bus/{sample}/transcripts.txt` | TXT | Transcript names |
| `index/{genome}.idx` | IDX | Kallisto index |

#### Salmon

The `salmon` directory contains output files generated by the `index` and `alevin` commands from the <ins>[Salmon](https://salmon.readthedocs.io/en/latest/)</ins> software:

| File | Format | Description |
| --- | --- | --- |
| `alevin/{sample}/quants_mat.gz` | TSV | Compressed count matrix |
| `alevin/{sample}/quants_mat_cols.txt` | TXT | Column header (gene_id) of the matrix |
| `alevin/{sample}/quants_mat_rows.txt` | TXT | Row index (CB-ids) of the matrix |
| `alevin/{sample}/quants_tier_mat.gz` | TSV | Tier categorization of the matrix |
| `index/{genome}` | DIR | Salmon index |

#### SingleCellExperiment

The `singlecellexperiment` directory contains <ins>[SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment/)</ins> object files containing spliced and unspliced count matrices:

| File | Format | Description |
| --- | --- | --- |
| `kallisto.rds` | RDS | SingleCellExperiment object from kallisto\|bustools workflow |
| `salmon.rds` | RDS | SingleCellExperiment object from Alevin workflow |
| `star.rds` | RDS | SingleCellExperiment object from STARsolo workflow |

#### STAR

The `star` directory contains output files generated by the `genomeGenerate` and `alignReads` commands from the <ins>[STAR](https://www.ncbi.nlm.nih.gov/pubmed/23104886)</ins> software:

| File | Format | Description |
| --- | --- | --- |
| `align/{sample}/Aligned.sortedByCoord.out.bam` | BAM | Spliced alignments |
| `align/{sample}/Solo.out/Gene` | DIR | Gene directory |
| `align/{sample}/Solo.out/Velocyto` | DIR | Velocyto directory |
| `index/{genome}` | DIR | STAR index |



## Tests

Test cases are in the `.test` directory. They are automatically executed via continuous integration with GitHub Actions.


## FAQ

#### Which quantification workflow should I use?

I would suggest there is no best workflow, each one captures a unique aspect of the data by the counting strategies they have implemented. For a more in-depth discussion, please refer to this research article by Soneson and colleagues: https://doi.org/10.1371/journal.pcbi.1008585

#### What reference genomes are supported?

The workflow uses genomepy and gffread to download and parse the user-specified reference genome and annotation. Therefore, any genome release compatible with these software should be supported.

#### How do I combine multiple sequencing runs?

If the sequencing runs were performed across multiple lanes on the same date, it is unlikely that a batch effect is present and I would recommend quantifying the files all together. Below is an example units table showing how to specify multiple sequencing runs jointly for a given sample: 

```
sample,unit,read1,read2
S1,L001,S1_L001.fastq.gz,S1_L001.fastq.gz
S1,L002,S1_L002.fastq.gz,S1_L002.fastq.gz
```

Alternatively, if the sequencing runs were performed on different machines and different dates, there is potential for a batch effect and I would recommend quantifying the files separately until this can be investigated. Below is an example units table showing how to specify multiple sequencing runs independently for a given sample:

```
sample,unit,read1,read2
S1_L001,L001,S1_L001.fastq.gz,S1_L001.fastq.gz
S1_L002,L002,S1_L002.fastq.gz,S1_L002.fastq.gz
```


## References

Chromium relies on multiple open-source software. Please give appropriate credit
by citing them in your publication:

[**Alevin**](https://doi.org/10.1186/s13059-019-1670-y)  
Srivastava, A., Malik, L., Smith, T. et al. Alevin efficiently estimates accurate gene abundances from dscRNA-seq data. Genome Biol 20, 65 (2019).

[**Anaconda**](https://anaconda.com)  
Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

[**BUSpaRse**](https://bioconductor.org/packages/BUSpaRse/)  
Moses L, Pachter L (2021). BUSpaRse: kallisto | bustools R utilities. R package version 1.4.2.

[**BUStools**](https://doi.org/10.1038/s41587-021-00870-2)  
Melsted, P., Booeshaghi, A.S., Liu, L. et al. Modular, efficient and constant-memory single-cell RNA-seq preprocessing. Nat Biotechnol (2021).

[**Bioconda**](https://doi.org/10.1038/s41592-018-0046-7)  
Grüning, B., Dale, R., Sjödin, A. et al. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods 15, 475–476 (2018).

[**eisaR**](https://bioconductor.org/packages/eisaR/)  
Stadler MB, Gaidatzis D, Burger L, Soneson C (2020). eisaR: Exon-Intron Split Anaalysis (EISA) in R. R package version 1.0.

[**genomepy**](http://dx.doi.org/10.21105/joss.00320)  
Heeringen, (2017), genomepy: download genomes the easy way, Journal of Open Source Software, 2(16), 320.

[**GffRead**](https://doi.org/10.12688/f1000research.23297.2)  
Pertea G and Pertea M. GFF Utilities: GffRead and GffCompare [version 2; peer review: 3 approved]. F1000Research 2020, 9:304.

[**Kallisto**](https://doi.org/10.1038/nbt.3519)  
Bray, N., Pimentel, H., Melsted, P. et al. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol 34, 525–527 (2016).

[**Python**](https://www.python.org/)  
Python Core Team (2015). Python: A dynamic, open source programming language.
Python Software Foundation.

[**R**](https://www.R-project.org/)  
R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.

[**STAR**](https://doi.org/10.1101/2021.05.05.442755)  
Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras, STAR: ultrafast universal RNA-seq aligner, Bioinformatics, Volume 29, Issue 1, January 2013, Pages 15–21.

[**SingleCellExperiment**](https://doi.org/10.1038/s41592-019-0654-x)  
Amezquita, R.A., Lun, A.T.L., Becht, E. et al. Orchestrating single-cell analysis with Bioconductor. Nat Methods 17, 137–145 (2020).

[**Snakemake**](https://doi.org/10.12688/f1000research.29032.2)  
Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 2; peer review: 2 approved]. F1000Research 2021, 10:33.
