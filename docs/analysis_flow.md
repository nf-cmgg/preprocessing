# nf-cmgg/preprocessing: Analysis Flow

## FGUMI-aware DNA flow (new branch)

The diagram below summarizes the current fgumi-aware branch and how it joins the common FASTQ to CRAM path.
It reflects the current wiring in `FASTQ_TO_CRAM` and `UMI_CONSENSUS_FGUMI`, including the early branching and data joins.

```mermaid
flowchart TD
  A[Input channel: meta + reads + aligner + index + fasta + gtf] --> B{sample_type}
  B -->|RNA| R1[FASTQ_ALIGN_RNA]
  B -->|DNA| C{meta.fgumi_aware == true}

  C -->|false| D1[FASTQ_ALIGN_DNA non-UMI]
  C -->|true| U0[Enter UMI_CONSENSUS_FGUMI]

    subgraph UMI_CONSENSUS_FGUMI
      U1[Step 1: FGUMI_EXTRACT\n(reads -> unmapped BAM with UMI tags)]
      U1J[Join with reference assets\nSNAP index + fasta + dict from meta.genome_data]
      U2[Step 3: FGUMI_SNAP_ZIPPER_SORT]
      U2a[samtools sort -n\nunmapped BAM]
      U2b[fgumi fastq]
      U2c[snap-aligner paired]
      U2d[samtools sort -n\npost-SNAP]
      U2e[fgumi zipper]
      U2f[fgumi sort --order template-coordinate]
      U3[Step 4: FGUMI_GROUP]
      U4[Step 5: FGUMI_SIMPLEX]
      U5J[Join simplex BAM with fasta]
      U5[Step 7a: FGUMI_FILTER]
      U6[Step 7b: FGUMI_SORT\ncoordinate sort + index]

      U1 --> U1J --> U2
      U2 --> U2a --> U2b --> U2c --> U2d --> U2e --> U2f
      U2f --> U3 --> U4 --> U5J --> U5 --> U6
    end

    D1 --> M1[Markdup branch selector\nbamsormadup | samtools | sort]
    R1 --> M1

    U6 --> MIX1[Mix UMI BAM/BAI into common postprocess stream]
    U3 --> MET1[grouping_metrics]
    U3 --> MET2[family_size_histogram]
    U4 --> MET3[consensus_metrics]
    U5 --> MET4[filtering_metrics]
    U6 --> MET5[filtered_consensus_bam]

    M1 --> P1[BIOBAMBAM_BAMSORMADUP or SAMTOOLS_SORMADUP or SAMTOOLS_SORT]
    P1 --> COMP{bam or cram}
    MIX1 --> COMP

    COMP -->|bam| CVT[SAMTOOLS_CONVERT to CRAM]
    COMP -->|cram| OUT1[cram + crai]
    CVT --> OUT1

    OUT1 --> E1[emit: cram_crai]
    MET1 --> E2[emit: sormadup_metrics]
    MET2 --> E3[emit: family_size_histogram]
    MET5 --> E4[emit: filtered_consensus_bam]
```

## Notes

- The fgumi branch is opt-in per sample via fgumi_aware.
  - Step numbering mirrors the implementation comments: steps 1, 3, 4, 5, and 7 are executed in this branch.
  - The "step 2" nomenclature from upstream fgumi Basic Workflow is intentionally absent in this pipeline path.
- UMI-specific metrics are emitted and mixed into the common reporting stream.
- The output still converges to the shared CRAM/crai downstream path.
