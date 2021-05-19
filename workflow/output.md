## Output <img align="right" width="200" src="../images/roundel.png">

### Contents

- [Overview](#overview)
- [Directory](#directory)
- [Files](#files)

### Overview

Chromium processes the data using three different quantification methods:

* kallisto|bustools
* Alevin
* STARsolo

For each method, all output files are saved to disk and a SingleCellExperiment  
object containing spliced and unspliced count matrices is generated.

### Directory

Chromium saves all output files in the `output` directory. The files are organised  
by software and named by sample wildcard and data type:

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

### Files

Chromium generates output files by software

#### BUSpaRse

The `busparse` directory contains output files generated by the `get_velocity_files` and `read_velocity_files` functions from the [BUSpaRse](https://bioconductor.org/packages/BUSpaRse/) software:

| File | Format | Description |
| --- | --- | --- |
| `{genome}.cDNA_introns.fa` | FASTA | Spliced transcript and intron sequences |
| `{genome}.cDNA_tx_to_capture.txt` | TXT | Transcript IDs of spliced transcripts |
| `{genome}.introns_tx_to_capture.txt` | TXT | Transcript IDs of introns |
| `{genome}.tr2g.tsv` | TSV | Transcripts and introns to genes |

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