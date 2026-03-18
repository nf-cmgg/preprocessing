# UMI consensus implementation (KAPA/fgbio)

This document describes the UMI consensus implementation used by the local `UMI_CONSENSUS_KAPA` workflow.

## Location

- Workflow: `subworkflows/local/umi_consensus/main.nf`
- Local process definitions: `modules/local/umi_consensus/main.nf`
- Integration point: `subworkflows/local/fastq_to_aligned_cram/main.nf`
- Process configuration: `conf/modules.config`

## Purpose

When `umi_consensus: true` is set on a non-RNA sample, the pipeline runs a UMI-aware consensus path before producing final sorted/indexed BAM output (later converted to CRAM in the parent subworkflow).

## Input and outputs

### Workflow input

`UMI_CONSENSUS_KAPA` expects channel entries shaped as:

- `[meta, reads, aligner, index, fasta]`

### Workflow outputs

- `bam_bai`: `[meta, bam, bai]`
- `family_sizes`: `[meta, histogram]`

## Process flow

1. `FASTQ_ALIGN_DNA`
   - Initial mapping of original reads with sample-selected aligner/index.
2. `UMI_LOCAL_SAMTOOLS_VIEW`
   - Local pre-filter step (`samtools view -F 260 -bh`) to keep primary mapped reads.
3. `UMI_SAMTOOLS_COLLATE` (nf-core module)
4. `UMI_SAMTOOLS_FIXMATE` (nf-core module)
5. `UMI_SAMTOOLS_SORT_TEMPLATE` (nf-core module)
   - Template-coordinate sort and index creation.
6. `UMI_FGBIO_COPYUMIFROMREADNAME` (nf-core module)
   - Copies UMI from read name to `RX` tag.
7. `UMI_FGBIO_GROUPREADSBYUMI` (nf-core module)
   - Groups reads into UMI families (`Adjacency` strategy).
8. `UMI_FGBIO_CALLMOLECULARCONSENSUSREADS` (nf-core module)
9. `UMI_FGBIO_FILTERCONSENSUSREADS` (nf-core module)
10. `UMI_SAMTOOLS_FASTQ` (nf-core module)
    - Converts filtered consensus BAM to interleaved FASTQ.
11. `FASTQ_ALIGN_DNA_CONSENSUS`
    - Re-maps consensus reads using original aligner/index.
12. `UMI_FGBIO_ZIPPERBAMS` (nf-core module)
    - Reconciles mapped consensus BAM with unmapped metadata/reference context.
13. `UMI_SAMTOOLS_SORT_FINAL` (nf-core module)
    - Final sort + index.

## Reference channel handling

The workflow derives helper reference channels from `meta`:

- `ch_meta_fasta_fai`: `[meta, fasta, fai]`
- `ch_meta_fasta`: `[meta, fasta]`
- `ch_meta_dict`: `[meta, dict]`

For sparse metadata in stubs/tests, `fai` and `dict` fallback to `/dev/null` to keep channel shapes valid.

## Config knobs

UMI process defaults are defined in `conf/modules.config` under selectors matching `.*FASTQ_TO_CRAM:UMI_.*` and module-specific names, including:

- `UMI_SAMTOOLS_COLLATE`
- `UMI_SAMTOOLS_FIXMATE`
- `UMI_SAMTOOLS_SORT_TEMPLATE`
- `UMI_FGBIO_GROUPREADSBYUMI`
- `UMI_FGBIO_CALLMOLECULARCONSENSUSREADS`
- `UMI_FGBIO_FILTERCONSENSUSREADS`
- `UMI_FGBIO_ZIPPERBAMS`
- `UMI_SAMTOOLS_SORT_FINAL`

## Why one local process remains

`UMI_LOCAL_SAMTOOLS_VIEW` is intentionally local. The upstream `samtools/view` module requires additional index/reference/qname inputs that do not fit this lightweight pre-filter use case without extra channel plumbing and collision risk.

## Testing

UMI process/workflow tests live under:

- `tests/modules/local/umi_consensus/`

The current stubs cover each process and the full workflow path.
