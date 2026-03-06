# nf-cmgg/preprocessing

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new/nf-cmgg/preprocessing)
[![GitHub Actions CI Status](https://github.com/nf-cmgg/preprocessing/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-cmgg/preprocessing/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-cmgg/preprocessing/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-cmgg/preprocessing/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-cmgg/preprocessing)

## Introduction

**nf-cmgg/preprocessing** is a bioinformatics pipeline that demultiplexes and aligns raw sequencing data.
It also performs basic QC and coverage analysis.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

Steps include:

- Demultiplexing using [`BCLconvert`](https://emea.support.illumina.com/sequencing/sequencing_software/bcl-convert.html)
- Run QC using [`MultiQC SAV`](https://github.com/MultiQC/MultiQC_SAV)
- Read QC and trimming using [`fastp`](https://github.com/OpenGene/fastp) or [`falco`](https://github.com/smithlabcode/falco)
- Alignment using either [`bwa`](https://github.com/lh3/bwa), [`bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2), [`bowtie2`](https://github.com/BenLangmead/bowtie2), [`dragmap`](https://github.com/Illumina/DRAGMAP), [`snap`](https://github.com/amplab/snap) or [`strobe`](https://github.com/ksahlin/strobealign) for DNA-seq and [`STAR`](https://github.com/alexdobin/STAR) for RNA-seq
- Duplicate marking using [`bamsormadup`](https://gitlab.com/german.tischler/biobambam2) or [`samtools markdup`](http://www.htslib.org/doc/samtools-markdup.html)
- Coverage analysis using [`mosdepth`](https://github.com/brentp/mosdepth) and [`samtools coverage`](http://www.htslib.org/doc/samtools-coverage.html)
- Alignment QC using [`samtools flagstat`](http://www.htslib.org/doc/samtools-flagstat.html), [`samtools stats`](http://www.htslib.org/doc/samtools-stats.html), [`samtools idxstats`](http://www.htslib.org/doc/samtools-idxstats.html) and [`picard CollectHsMetrics`](https://broadinstitute.github.io/picard/command-line-overview.html#CollectHsMetrics), [`picard CollectWgsMetrics`](https://broadinstitute.github.io/picard/command-line-overview.html#CollectWgsMetrics), [`picard CollectMultipleMetrics`](https://broadinstitute.github.io/picard/command-line-overview.html#CollectMultipleMetrics)
- QC aggregation using [`multiqc`](https://multiqc.info/)

<picture>

  <source media="(prefers-color-scheme: dark)" srcset="docs/images/metro_map_dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="docs/images/metro_map_light.svg">
  <img alt="Fallback image description" src="docs/images/metro_map_light.svg">
</picture>

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

The full documentation can be found [here](docs/README.md)

First, prepare a samplesheet with your input data. Check the [usage docs](docs/usage.md) for details on the required format and example files.

Now, you can run the pipeline using:

```bash
nextflow run nf-cmgg/preprocessing \
   -profile <docker/singularity/...> \
   --igenomes_base /path/to/genomes \
   --input samplesheet.<csv|yaml|json> \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

nf-cmgg/preprocessing was originally written by the CMGG ICT team.

## Support

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
