# CenterForMedicalGeneticsGhent/nf-cmgg-preprocessing: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0dev

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
