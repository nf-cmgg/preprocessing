# ![nf-core/cmggpreprocessing](docs/images/nf-core/cmggpreprocessing_logo_light.png#gh-light-mode-only) ![nf-core/cmggpreprocessing](docs/images/nf-core/cmggpreprocessing_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/cmggpreprocessing/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/cmggpreprocessing/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/cmggpreprocessing/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/cmggpreprocessing/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/cmggpreprocessing/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23cmggpreprocessing-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/cmggpreprocessing)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**nf-core/cmggpreprocessing** is a bioinformatics best-practice analysis pipeline for Preprocessing workflow for sequencing data at CMGG.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/cmggpreprocessing/results).

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run nf-core/cmggpreprocessing -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```console
   nextflow run nf-core/cmggpreprocessing --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Documentation

The nf-core/cmggpreprocessing pipeline comes with documentation about the pipeline [usage](https://nf-co.re/cmggpreprocessing/usage), [parameters](https://nf-co.re/cmggpreprocessing/parameters) and [output](https://nf-co.re/cmggpreprocessing/output).

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
    DEMUX_FASTQ                                 --> FASTP[FastP]
    FASTP                                       --> DEMUX_MULTIQC[MultiQC]
    DEMULTIPLEX_STATS                           --> DEMUX_MULTIQC[MultiQC]
    DEMUX_MULTIQC                               --> DEMUX_MULTIQC_REPORT([MultiQC Report])
end

DEMULTIPLEX                                     --> RAW_FASTQ([Demultiplexed Fastq per sample per lane])
DEMULTIPLEX                                     --> DEMUX_REPORTS([Demultiplexing reports])
RAW_FASTQ                                       --> ALIGNMENT

subgraph ALIGNMENT
    direction TB
    FASTQ([Fastq per sample per lane])          --> IS_HUMAN{Human data?}
    IS_HUMAN                                    --YES--> NX_SPLITFASTQ[nx SplitFastq]
    IS_HUMAN                                    --NO--> FASTQTOSAM[Picard FastqToSam]
    FASTQTOSAM                                  --> UNALIGNED_BAM([Unaligned BAM])
    NX_SPLITFASTQ                               --split by # reads--> ALIGNER{ALIGNER}
    ALIGNER                                     --> BOWTIE2[Bowtie2-align]
    ALIGNER                                     --> BWAMEM2[Bwamem2 mem]
    ALIGNER                                     --> SNAP[snap-aligner]
    BOWTIE2                                     --> MERGE((Merge bam chunks))
    BWAMEM2                                     --> MERGE((Merge bam chunks))
    SNAP                                        --> MERGE((Merge bam chunks))
    MERGE                                       --> BAMPROCESSOR{BAMPROCESSOR}
    BAMPROCESSOR{BAMPROCESSOR}                  --> BIOBAMBAM[Biobambam]
    BAMPROCESSOR{BAMPROCESSOR}                  --> ELPREP[Elprep]

    subgraph ELPREP_FLOW[Elprep subroutine]
        direction TB
        ELPREP_BAM([BAM])                       --> ELPREP_SPLIT[Elprep split]
        ELPREP_SPLIT                            --split by chromosome-->    ELPREP_FILTER[Elprep filter]
        ELPREP_FILTER                           --sort/mark duplicates-->   ELPREP_MERGE[Elprep merge]
        ELPREP_FILTER                           --BQSR/variant calling-->   ELPREP_MERGE[Elprep merge]
        ELPREP_MERGE                            --> ELPREP_SORTBAM([Postprocessed BAM])
        ELPREP_MERGE                            --> ELPREP_GVCF([gVCF])
    end
    ELPREP                                      --> ELPREP_FLOW
    ELPREP_FLOW                                 --> SORTBAM([Postprocessed BAM])
    ELPREP_FLOW                                 --> ELPREP_GVCF([gVCF])
    ELPREP_FLOW                                 --> MARKDUP_METRICS([Markduplicates Metrics])

    subgraph BIOBAMBAM_FLOW[Biobambam subroutine]
        direction TB
        BIOBAMBAM_BAM([BAM])                    --> SPLIT[Split tool TBD]
        SPLIT                                   --split by chromosome-->    BAMSORMADUP[BamSorMaDup]
        BAMSORMADUP                             --sort/mark duplicates -->  BIOBAMBAM_SORTBAM([Postprocessed BAM])
    end
    BIOBAMBAM                                   --> BIOBAMBAM_FLOW
    BIOBAMBAM_FLOW                              --> SORTBAM
    BIOBAMBAM_FLOW                              --> MARKDUP_METRICS([Markduplicates Metrics])

    SORTBAM                                     -->  BAMQC[BAM QC Tools]
    BAMQC                                       -->  BAM_METRICS([BAM metrics])
    SORTBAM                                     -->  SCRAMBLE[Scramble]
    SORTBAM                                     -->  MOSDEPTH[Mosdepth]
    MOSDEPTH                                    -->  COVERAGE_BED([Coverage BEDs])
    MOSDEPTH                                    -->  COVERAGE_METRICS([Coverage Metrics])
    UNALIGNED_BAM                               -->  SCRAMBLE
    SCRAMBLE                                    -->  ALN_CRAM([CRAM])

    BAM_METRICS                                 -->  ALN_MULTIQC[MultiQC]
    MARKDUP_METRICS                             -->  ALN_MULTIQC
    COVERAGE_METRICS                            -->  ALN_MULTIQC
    ALN_MULTIQC                                 -->  ALN_MULTIQC_REPORT([MultiQC Report])
end

ALIGNMENT                                       --> A_CRAM([CRAM])
ALIGNMENT                                       --> GVCF([gVCF])
ALIGNMENT                                       --> A_BAM_METRICS([BAM metrics])
ALIGNMENT                                       --> A_COVERAGE_METRICS([Coverage metrics])
ALIGNMENT                                       --> A_COVERAGE_BED([Coverage BED])

A_BAM_METRICS                                   --> MQC([Run MultiQC Report])
A_COVERAGE_METRICS                              --> MQC
DEMUX_REPORTS                                   --> MQC

```


## Credits

nf-core/cmggpreprocessing was originally written by Matthias De Smet.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#cmggpreprocessing` channel](https://nfcore.slack.com/channels/cmggpreprocessing) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/cmggpreprocessing for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
