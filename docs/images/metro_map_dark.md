```mermaid
%%metro logo: ./nf-cmgg-preprocessing_logo_dark.png
%%metro style: dark
%%metro line: main | Alignment and Postprocessing | #00ff00
%%metro line: qc | Quality control | #ff0000
%%metro file: BCL_IN | BCL
%%metro file: FASTQ_IN | FASTQ
%%metro file: CRAM_OUT | CRAM
%%metro file: MULTIQC_LIBRARY | HTML
%%metro file: MULTIQC_SAV | HTML
%%metro compact_offsets: true

graph TD

    BCL_IN[]
    FASTQ_IN[]
    MULTIQC_SAV[]
    CRAM_OUT[]
    MULTIQC_LIBRARY[]

    BCL_IN -->|main | BCLCONVERT
    BCL_IN -->|qc| MULTIQC_SAV
    BCLCONVERT -->|qc| MULTIQC_SAV

    FASTQ_IN[]
    FASTQ_IN -->|qc,main| FASTP
    BCLCONVERT -->|qc| FALCO
    BCLCONVERT -->|qc,main| FASTP
    FALCO -->|qc| MULTIQC_LIBRARY
    FASTP -->|qc| MULTIQC_LIBRARY

    FASTP -->|main| UMI_CONSENSUS
    FASTP -->|main| ALIGN
    UMI_CONSENSUS -->|main| ALIGN
    ALIGN -->|main| MARKDUP
    MARKDUP -->|main| CRAM_OUT

    CRAM_OUT -->|qc| MOSDEPTH
    CRAM_OUT -->|qc| SAMTOOLS_COV
    MOSDEPTH -->|qc| MULTIQC_LIBRARY
    SAMTOOLS_COV -->|qc| MULTIQC_LIBRARY

    CRAM_OUT -->|qc| SAMTOOLS_QC
    CRAM_OUT -->|qc| PICARD
    SAMTOOLS_QC -->|qc| MULTIQC_LIBRARY
    PICARD -->|qc| MULTIQC_LIBRARY
```
