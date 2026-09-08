# nf-cmgg/preprocessing: Usage

Parameter documentation can be found [here](parameters.md)

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It can be a CSV, TSV, JSON or YAML file.

```bash
--input '[path to samplesheet file]'
```

The pipeline supports two types of samplesheets to be used as input: [`fastq`](#fastq-samplesheet) and [`flowcell`](#flowcell-samplesheet) samplesheets. The type will be automatically detected and applied by the pipeline. FASTQ rows need `fastq_1`; `fastq_2` is optional for single-end data. After demultiplexing, single-end FASTQs are handled if only one read file is produced.

### Fastq samplesheet

A `fastq` samplesheet file consisting of paired-end data may look something like the one below.

```yml
- id: DNA1_L001
  samplename: DNA_paired1
  library: test_library
  genome: GRCh38
  aligner: bwamem
  markdup: bamsormadup
  umi_aware: false
  skip_trimming: false
  trim_front: 0
  trim_tail: 0
  adapter_R1: AGATCGGAAGAGCACACGTCTGAACTCCTTA
  adapter_R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
  qc_mode: basic
  roi: null
  tag: WES
  sample_type: DNA
  fastq_1: https://github.com/nf-cmgg/test-datasets/raw/preprocessing/data/genomics/homo_sapiens/illumina/fastq/sample1_R1.fastq.gz
  fastq_2: https://github.com/nf-cmgg/test-datasets/raw/preprocessing/data/genomics/homo_sapiens/illumina/fastq/sample1_R2.fastq.gz
```

Following table shows the fields that are used by the `fastq` samplesheet:

| Column                               | Description                                                                                                                                                                                                                        | Required                                        |
| ------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------- |
| `id`                                 | Unique sample identifier                                                                                                                                                                                                           | :heavy_check_mark:                              |
| `samplename`                         | The sample name corresponding to the sample in the Fastq file(s)                                                                                                                                                                   | :heavy_check_mark:                              |
| `genome`                             | Genome build. Allowed values: `GRCh38`, `GRCh38-noalt`, `GRCm39`, `GRCz11`, `hg38`, `hg38-noalt`. See [organism to genome mapping](#organism-to-genome-mapping) for what is used when only `organism` is set.                      | :heavy_check_mark: (unless `organism` is given) |
| `organism`                           | Full name of the organism. Currently supports `Homo sapiens`, `Mus musculus`, `Danio rerio` and `Equus caballus`                                                                                                                   | :heavy_check_mark: (unless `genome` is given)   |
| `library`                            | Sample library name. When set, results are published under `library/samplename`.                                                                                                                                                   | :x:                                             |
| `tag`                                | Sample tag (`[A-Za-z0-9_-]+`). `SeqCap` restricts panel coverage gene lists to files whose names contain `seqcap`.                                                                                                                 | :x:                                             |
| `aligner`                            | DNA aligner: `bowtie2`, `bwamem`, `bwamem2`, `dragmap`, `strobe` or `snap`. Set to `false` to skip alignment and emit FASTQ. RNA samples (`sample_type: RNA`) always use `star`.                                                   | :heavy_check_mark:                              |
| `markdup`                            | Markdup algorithm to use for duplicate marking. Can be set to `bamsormadup`, `samtools` or `false`                                                                                                                                 | :x:                                             |
| `umi_aware`                          | Whether UMI-aware processing should be used. Only applies when `markdup` is set to `samtools`                                                                                                                                      | :x:                                             |
| `call_consensus`                     | Perform consensus calling using the `fgumi` toolsuite. This only works for DNA samples and will always run the SNAP aligner                                                                                                        | :x:                                             |
| `fgumi_simplex_min_reads`            | Minimum number of reads required per UMI family for fgumi simplex consensus generation. Defaults to `1` and should be `1` or higher.                                                                                               | :x:                                             |
| `fgumi_snap_ignore_mismatched_pairs` | Pass -I to SNAP to ignore mismatched read IDs in paired-end input when using the `fgumi` toolsuite (`call_consensus` set as `true`). Defaults to `true`                                                                            | :x:                                             |
| `skip_trimming`                      | Skip adapter trimming step                                                                                                                                                                                                         | :x:                                             |
| `trim_front`                         | Number of bases to trim from the front of reads                                                                                                                                                                                    | :x:                                             |
| `trim_tail`                          | Number of bases to trim from the tail of reads                                                                                                                                                                                     | :x:                                             |
| `adapter_R1`                         | Adapter sequence for read 1                                                                                                                                                                                                        | :x:                                             |
| `adapter_R2`                         | Adapter sequence for read 2                                                                                                                                                                                                        | :x:                                             |
| `qc_mode`                            | QC mode for the sample. Can be set to `basic` or `full`. Basic QC includes samtools flagstat, idxstats and mosdepth. Full QC includes samtools stats, samtools coverage, riker metrics and panel coverage in addition to basic QC. | :x:                                             |
| `roi`                                | The path to a BED file containing Regions Of Interest for coverage analysis                                                                                                                                                        | :x:                                             |
| `sample_type`                        | Sample type. Allowed values: `DNA`, `RNA`, `Tissue`. Defaults to `DNA`. RNA samples are aligned with STAR.                                                                                                                         | :x:                                             |
| `fastq_1`                            | FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'                                                                                                                     | :heavy_check_mark:                              |
| `fastq_2`                            | FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'. Omit for single-end data.                                                                                                            | :x:                                             |

An [example samplesheet](../tests/inputs/test.yml) has been provided with the pipeline.

### Flowcell samplesheet

A `flowcell` samplesheet file consisting of one sequencing run may look something like the one below.

```yml
- id: 200624_A00834_0183_BHMTFYDRXX
  samplesheet: https://github.com/nf-cmgg/test-datasets/raw/refs/heads/preprocessing/data/genomics/homo_sapiens/illumina/flowcell/SampleSheet_2.csv
  lane: 1
  flowcell: s3://test-data/genomics/homo_sapiens/illumina/bcl/
  sample_info: https://github.com/nf-cmgg/test-datasets/raw/refs/heads/preprocessing/data/genomics/homo_sapiens/illumina/flowcell/SampleInfo_2.json
```

Following table shows the fields that are used by the `flowcell` samplesheet:

| Column        | Description                                                                                                  | Required           |
| ------------- | ------------------------------------------------------------------------------------------------------------ | ------------------ |
| `id`          | Unique flowcell identifier                                                                                   | :heavy_check_mark: |
| `samplesheet` | Illumina sample sheet CSV for the flowcell lane                                                              | :heavy_check_mark: |
| `sample_info` | JSON/YAML file with sample information. See the [flowcell sample info](#flowcell-sample-info) documentation. | :heavy_check_mark: |
| `flowcell`    | Illumina flowcell directory                                                                                  | :heavy_check_mark: |
| `lane`        | Lane number                                                                                                  | :x:                |

An [example samplesheet](../tests/inputs/test.yml) has been provided with the pipeline.

### Flowcell sample info

A `flowcell` sample info JSON/YAML file for one sequencing run may look something like the one below.

```yml
- samplename: DNA_paired1
  library: test_library
  genome: GRCh38
  aligner: bwamem
  markdup: bamsormadup
  umi_aware: false
  skip_trimming: false
  trim_front: 0
  trim_tail: 0
  adapter_R1: AGATCGGAAGAGCACACGTCTGAACTCCTTA
  adapter_R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
  qc_mode: basic
  roi: null
  tag: WES
  sample_type: DNA
```

Each row needs `samplename`, `aligner`, `tag`, and either `genome` or `organism`. Analysis fields match the [fastq samplesheet](#fastq-samplesheet). Extra sample-info fields are `purpose` (`research` or `diagnostic`), `vivar_project`, `binsize`, and `panels`.

The same `samplename` may appear in more than one library. In that case each row needs a distinct `library` value, and the Illumina sample sheet must set `LibraryName` so that BCL Convert `RGLB` matches `sampleinfo.library`. A single row per `samplename` does not need `RGLB` on the demultiplexed FASTQ.

### Organism to genome mapping

When a row sets `organism` but no `genome`, the pipeline derives the genome build:

| Organism         | Genome build |
| ---------------- | ------------ |
| `Homo sapiens`   | `GRCh38`     |
| `Mus musculus`   | `mm10`       |
| `Danio rerio`    | `GRCz11`     |
| `Equus caballus` | `EquCab2`    |

Matching is case-insensitive and also accepts an underscore instead of a space. Any other organism leaves the genome unset. Samples whose genome has no entry in `conf/igenomes.config` skip alignment and are QC'd with falco.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-cmgg/preprocessing --input ./samplesheet.<csv|json|yaml> --outdir ./results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-cmgg/preprocessing -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-cmgg/preprocessing
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-cmgg/preprocessing releases page](https://github.com/nf-cmgg/preprocessing/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Apple containers) - see below.

The pipeline loads institutional configs from [nf-core/configs](https://github.com/nf-core/configs) by default (`params.custom_config_base`). If `custom_config_base` points at an nf-cmgg configs tree, it also loads `pipeline/preprocessing.config` from that repo.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `debug`
  - A generic profile with settings to help with debugging the pipeline. It will use more verbose logging.
- `arm64`
  - A generic profile with settings to run the pipeline on ARM64 architecture machines (eg. Apple Silicon). It will use software containers built for ARM64 where available.
- `emulate_amd64`
  - A generic profile with settings to run the pipeline on ARM64 architecture machines (eg. Apple Silicon) using AMD64 software containers. This is for when ARM64 containers are not available but you still want to run the pipeline on an ARM64 machine. Note that this will be slower than using ARM64 containers.
- `apple`
  - A generic configuration profile to be used with [Apple containers](https://github.com/apple/container)
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `test_full`
  - A profile with a more complete test dataset
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - Enable Seqera Wave to resolve containers
- `gpu`
  - Pass GPU flags through to Docker, Apptainer or Singularity

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version may be out of date.

To use a different container from the default container specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### Institutional configs

Override `custom_config_base` (and usually `custom_config_version`) to use [nf-cmgg/configs](https://github.com/nf-cmgg/configs) instead of the default nf-core configs. Test a one-off config with `-c` first.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted to your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~/.bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
