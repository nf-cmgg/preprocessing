# CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing

[![GitHub Actions CI Status](https://github.com/CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing/workflows/nf-core%20CI/badge.svg)](https://github.com/CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing/workflows/nf-core%20linting/badge.svg)](https://github.com/CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing/actions?query=workflow%3A%22nf-core+linting%22)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)

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

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

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
   nextflow run CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing --input flowcells.csv --samples samples.csv --outdir <OUTDIR> -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Flowchart

```mermaid

flowchart TB

FC(["Flowcell (BCL)"])                          --> DEMULTIPLEX
SS([SampleSheet])                               --> DEMULTIPLEX

subgraph DEMULTIPLEX[Demultiplex]
    direction LR
    SAMPLESHEET([SampleSheet])                  --> BCLCONVERT[bcl-convert]
    FLOWCELL([Flowcell])                        --Split by LANE--> BCLCONVERT[bcl-convert]
    BCLCONVERT                                  --> DEMUX_FASTQ([Fastq])
    BCLCONVERT                                  --> DEMULTIPLEX_STATS([Demultiplex Reports])
end

DEMULTIPLEX                                     --> FASTP[FastP: Trimming and QC]
DEMULTIPLEX                                     --> DEMUX_REPORTS
FASTP                                           --> IS_HUMAN{Human data?}
IS_HUMAN                                        --YES--> ALIGNMENT
IS_HUMAN                                        --NO--> FASTQTOSAM
FASTQTOSAM[Picard FastqToSam]                   --> UNALIGNED_BAM([Unaligned BAM]) --> A_CRAM

subgraph ALIGNMENT
    direction TB

    subgraph ALIGNER
        direction LR
        BOWTIE2[bowtie2-align] & BWAMEM2[bwamem2 mem] & SNAP[snap-aligner] & DRAGMAP[dragmap] --> SORT[Sorting]
    end

    ALIGNER --> BamSorMaDUP
    BamSorMaDUP                                 --> SORTBAM[Sorted/markdup bam] & MARKDUP_METRICS([Markduplicates Metrics])
    SORTBAM                                     -->  ALN_CRAM([CRAM]) & MOSDEPTH[Mosdepth] & BAMQC[BAM QC Tools]
    MOSDEPTH                                    -->  COVERAGE_BED([Coverage BEDs]) & COVERAGE_METRICS([Coverage Metrics])
    BAMQC & MARKDUP_METRICS & COVERAGE_METRICS  -->  ALN_MULTIQC[MultiQC]
end

ALIGNMENT                                       --> A_CRAM([CRAM])
ALIGNMENT                                       --> A_BAM_METRICS([BAM metrics])
ALIGNMENT                                       --> A_COVERAGE_METRICS([Coverage metrics])
ALIGNMENT                                       --> A_COVERAGE_BED([Coverage BED])

A_BAM_METRICS                                   --> MQC([Run MultiQC Report])
A_COVERAGE_METRICS                              --> MQC
FASTP                                           --> MQC
DEMUX_REPORTS                                   --> MQC

```

## Credits

CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing was originally written by Matthias De Smet.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
