# CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing

[![GitHub Actions CI Status](https://github.com/CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing/workflows/nf-core%20CI/badge.svg)](https://github.com/CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing/workflows/nf-core%20linting/badge.svg)](https://github.com/CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing/actions?query=workflow%3A%22nf-core+linting%22)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)

## Introduction

**CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing** is a bioinformatics best-practice analysis pipeline for sequencing data preprocessing at CMGG.
This workflow handles demultiplexing, alignment, qc and archiving.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

1. Demultiplexing
2. Fastq preprocessing
3. Alignment
4. QC
5. Coverage analysis
6. Compression

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=23.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```console
   nextflow run CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing --input samplesheet<.csv/.yaml> --outdir <OUTDIR> -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Flowchart

```mermaid

flowchart TB

subgraph DEMULTIPLEX[Basecalling & Demultiplex]
    direction LR
    SampleSheet --> BCL-convert
    FlowCell    --Split by LANE--> BCL-convert
    BCL-convert --> DEMULTIPLEX_STATS([Demultiplex Reports])
    BCL-convert --> DEMUX_FASTQ([FastQ])
end

subgraph DEALIGNMENT[Bam/Cram to FastQ Conversion]
  direction TB
  BAM/CRAM --> SAMTOOLS_COLLATE[Samtools Collate] --> SAMTOOLS_FASTQ[Samtools Fastq] --> DEALIGNED_FASTQ([FastQ])
end

subgraph ALIGNMENT[Alignment]
    direction LR
    subgraph ALIGNER
      direction TB
      BWA-mem --> BAM
      BWAmem2 --> BAM
      DragMap --> BAM
      Bowtie2 --> BAM
      Snap    --> BAM
    end
    FQ_TO_ALIGN([FastQ]) --> ALIGNER
    ALIGNER --> Merge/Sort/MarkDuplicates --> Cram([Cram])
end

subgraph FQtoUCRAM[FastQ to Unaligned CRAM conversion]
    direction TB
    FQ_TO_CONVERT([FastQ])
    --> Samtools_import
    --> Unaligned_CRAM([Unaligned CRAM])
end

subgraph QC[Quality Assessment]
    direction TB
    CRAM_TO_QC([Cram]) --> SAMTOOLS --> stats & idxstats & Flagstat
    CRAM_TO_QC --> PICARD --> CollectWGSmetrics/CollectHsMetrics & CollectMultipleMetrics
end

subgraph COVERAGE[Coverage Analysis]
    direction TB
    CRAM_TO_COVERAGE([Cram]) --> Mosdepth --> BED([Coverage BED])
end

FQ_INPUT([FastQ Input]) --> SUPPORTED{Supported Genome?}

FC_INPUT([Flowcell Input]) --> DEMULTIPLEX --> SUPPORTED{Supported Genome?}
BAMCRAM_INPUT([Bam/Cram Input]) --> DEALIGNMENT --> SUPPORTED{Supported Genome?}

SUPPORTED{Supported Genome?} --> |Yes| ALIGNMENT
SUPPORTED{Supported Genome?} --> |No| FQtoUCRAM

ALIGNMENT --> QC
ALIGNMENT --> COVERAGE
COVERAGE --> MQC[MultiQC]
QC --> MQC

```

## Credits

CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing was originally written by Matthias De Smet.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
