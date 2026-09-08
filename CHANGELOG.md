# nf-cmgg/preprocessing: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 3.1.0 (dev)

- Fail DNA alignment when the requested aligner is not supported instead of dropping the sample.
- Keep UMI consensus filter/SNAP steps on the merged sample meta so they still run after fastp splitting.
- Omit empty FASTQ `CN`/`LB` read-group fields so they are not written as empty strings.
- Add support for Equus caballus genome (EquCab2) when the organism is specified in the sample metadata.
- Associate demultiplexed FASTQs with multiplexed `sampleinfo` rows using `readgroup.LB` when the same `samplename` appears in more than one library ([#169](https://github.com/nf-cmgg/preprocessing/issues/169)). Those libraries are published under `LIBRARY/SAMPLENAME`.
- Drop `Picard` modules and associated parameters in favor of `riker multi`, which is much more efficient and now run by default.
- Add `apple` profile to enable the use of Apple containers.
- Drop `hook_url` parameter. Notifications are now handled by the `nf-teams` plugin, which is enabled in `nextflow.config` and configured through a custom config.
- Drop `run_coverage` and add `qc_mode` parameter. Thorough coverage analysis is now included in the `qc_mode` parameter, which can be set to `basic` or `full`. Basic QC includes samtools flagstat, idxstats and mosdepth. Full QC includes samtools stats, samtools coverage, riker metrics and panel coverage in addition to basic QC.
- Revert default custom config to `nf-core` source
- Bump pipeline to nf-core template v4.1.0
- Bump modules to latest versions

## 3.0.2

- remove common plugins in favor of defining them in the nf-cmgg/configs, which will be used across all nf-cmgg pipelines. This allows for better version control and consistency across pipelines, as well as reducing the maintenance burden of keeping plugins up to date in multiple repositories.

## 3.0.1

- Add parameter typing to make strict syntax runs work more properly
- Bumped minimal Nextflow version to 26.04.0 to allow for strict syntax and other improvements
- Fix an issue with `mosdepth` environment variables not being correctly set on GCP
- Bump `mosdepth` module to 0.3.14
- Bump nf-core template to v4.0.2
- Fix fastq outputs when aligner is set to `false` or unsupported genome, which were not being correctly emitted in the previous version
- Fix `FALCO` module not being correctly run on fastq samples, which was caused by an issue with the branching logic in the previous version

## 3.0.0

- Add `MultiQC-SAV` module for Illumina Run QC reports
- Migrate readgroup parsing from `preprocessing.nf` to a local subworkflow
- Update pipeline to use topic channels only, deprecation `versions` channels
- Drop `bcl_demultiplex` subworkflow in favour of `bcl-convert` module
- Update the output handling to use the new workflow output definitions.
- Bump all modules to latest versions.
- The workflow now outputs data in a subdirectory per `library`, including a library specific MultiQC report
- Drop support for unaligned cram outputs in favor of untrimmed fastq outputs, which are more widely supported and can be used for a wider range of downstream analyses.
- Add support for untrimmed fastq outputs for unsupported genomes or when aligner is set to `false`.
- Drop support for global `aligner` parameter. The aligner must now be specified per sample in the sample sheet or sample info.
- Drop support for global `markdup` and `umi_aware` parameters. Marking duplicates must now be specified per sample in the sample sheet or sample info.
- Drop support for global `run_coverage` and `disable_picard_metrics` parameters. Running coverage analysis must now be specified per sample in the sample sheet or sample info.
- Drop support for global `skip_trimming`, `trim_front`, `trim_tail`, `adapter_R1` and `adapter_R2` parameters. Trimming must now be specified per sample in the sample sheet or sample info.
- Drop support for global `roi` parameter. Regions of interest must now be specified per sample in the sample sheet or sample info.
- Simplify fastq sharding and make it user configurable via the `split_fastq` parameter.
- Added splice junctions and junctions outputs for RNA-seq alignments using STAR.

## v2.0.6

- Update bash options for nextflow version 25.04.2
- Update `bcl-convert` module to 4.3.13
- Update `bowtie2` module to 2.5.4
- Update `multiqc` module to 1.29

## v2.0.5

- Panel coverage will now run all genelists' coverage analysis in one process per sample
- Fix readgroup sample names for `snp_` prefixed samples
- Edit `WGS` profile to include `samtools/sormadup` and enable UMI aware duplicate marking

## v2.0.4

- Add flag to `snap` to allow for more expansion of compressed data
- Parse readgroup samplenames for `snp_` prefixed samples
- Replace GRCh39 by mm10 as default genome for `Mus musculus`

## v2.0.3

- Update modules
- Fix bugs where samples with unequal number of files were not processed correctly
- Update to nf-core template v3.1.1
- Fix all tests

## v2.0.2

- Drop unused params
- Set aligner to `star` for RNA-seq
- Finetune resources
- Fix some bugs for different input tags
- Update modules

## v2.0.1

- Fix syntax according to new rules
- Drop usage of `check_max` function in favour of `resourceLimits`
- Bump nf-schema to v2.1.0
- Add configuration profiles for several datatypes
- Fix usage of `--run_coverage` parameter
- Update modules

## v2.0.0

- Move repo to nf-cmgg organisation
- Overhaul reference data handling. All data is now sourced from the `igenomes` config file and dynamically loaded based on the `organism` metadata field.
- Drop support for bam/cram inputs
- Update to nf-core template v2.13.1
- Drop `fgbio/fastqtosam` in favour of `samtools/import`
- Add `samtools/coverage` for coverage analysis
- Add testing with nf-test
- Overhaul panel coverage analysis
- Add support for `STAR` aligner for RNA-seq
- Replace CSV sampleinfo definition with validated JSON/YAML
- Update modules

## v1.2.0

- Add fix for sample replicates with different split sizes
- Add option to use samtools for duplicate marking
- Updated modules
- Updated to nf-core template v2.8
- Fix issue where unaligned samples were not merged and fixed naming
- Add minimal support for sample "tags" in sample sheet
- Add panel coverage analysis
- Limit panel coverage analyses to samples with WES/WGS tags
- Make BAM_ARCHIVE accept channels with cram files
- Drop ELPREP support
- Add support for cram and index outputs in SAMTOOLS_SORMADUP
- Add option to disable marking duplicates
- Update to nf-core template v2.9
- Start using nf-validation for input parsing
- Drop support for generating aligner indices

## v1.1.2

- Support for bamprocessing with ELPREP
- Improved sample grouping
- Meta value clean up
- bump modules
- fix multiqc report to include all results
- fix duplicate marking with bamsormadup
- drop post alignment sorting, doesn't make sense as we do sort/markdup afterwards
- Add option to use elprep for alignment postprocessing

## v1.0.2dev

- Add `unaligned` denominator to unaligned samples

## v1.0.1

- Updated modules

### `Fixed`

- Bug with `fgbio/fastqtosam` where having multiple files per read requires a read structure

## v1.0dev - [06/09/2022]

Initial release of CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing, created with the [nf-core](https://nf-co.re/) template.
