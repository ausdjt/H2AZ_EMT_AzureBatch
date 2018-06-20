# H2AZ_EMT
Bash and R scripts associated with analysis of H2A.Z ChIP-Seq data from TGFb-induced EMT in MDCK cells. Pre-release.

# some requirements and status of different pipeline stages
## tools & binaries
* snakemake [>4.0]
* bowtie2 [>2.1]
* cutadapt
* samtools
* deepTools
* R/Bioconductor [>3.4]
* picardTools
* htseq


## pipeline stages
* cutadapt_subworkflow.py - OK
* bowtie2_subworkflow.py - OK
* bam_processing_subworkflow.py - OK
* merge_replicates_subworkflow.py - seems OK?
* deepTools_QC_subworkflow.py - OK

## Workflow invocation
Following variables need to be set at the command line at time of invocation of the workflow:
* ASSAY="<type of sequencing experiment, e.g. ChIP-Seq or RNA-Seq, relative to working dir>"
* RUNID="<name of subdirectory containing sequencing data, e.g. ArrayExpress ID, relative to working dir>"
* WORKFLOWDIR="<name of directory where snakemake files are stores, e.g. name of git repo, relative to user /home>"

### Example snakemake command
```
snakemake --snakefile merge_replicates_subworkflow.py --configfile config.json --config ASSAY=ChIP-Seq RUNID=NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq WORKFLOWDIR="/Development/JCSMR-Tremethick-Lab/" --jobs 32 -pr
```
## Example directory structure
variables used:
ASSAY="ChIP-Seq"
RUNID="NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"
```
ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/
├── fastq
└── processed_data
    ├── CanFam3.1_ensembl84_ERCC
    │   ├── bowtie2
    │   │   ├── duplicates_marked
    │   │   ├── duplicates_removed
    │   │   ├── quality_filtered
    │   │   └── sorted
    │   ├── deepTools
    │   │   ├── bamCompare
    │   │   │   └── normal
    │   │   │       ├── duplicates_marked
    │   │   │       │   ├── readCount
    │   │   │       │   └── SES
    │   │   │       └── duplicates_removed
    │   │   │           ├── readCount
    │   │   │           └── SES
    │   │   ├── bamCoverage
    │   │   │   ├── debugging
    │   │   │   ├── MNase
    │   │   │   │   ├── duplicates_marked
    │   │   │   │   └── duplicates_removed
    │   │   │   └── normal
    │   │   │       ├── duplicates_marked
    │   │   │       └── duplicates_removed
    │   │   ├── bamPEFragmentSize
    │   │   │   ├── duplicates_marked
    │   │   │   └── duplicates_removed
    │   │   ├── bigwigCompare
    │   │   │   └── normal
    │   │   │       ├── duplicates_marked
    │   │   │       │   ├── readCount
    │   │   │       │   └── SES
    │   │   │       └── duplicates_removed
    │   │   │           ├── readCount
    │   │   │           └── SES
    │   │   ├── computeMatrix
    │   │   │   ├── reference-point
    │   │   │   │   ├── bigwigCompare
    │   │   │   │   │   ├── duplicates_marked
    │   │   │   │   │   │   └── TSS
    │   │   │   │   │   └── duplicates_removed
    │   │   │   │   │       └── TSS
    │   │   │   │   ├── duplicates_marked
    │   │   │   │   │   └── TSS
    │   │   │   │   └── duplicates_removed
    │   │   │   │       └── TSS
    │   │   │   └── scale-regions
    │   │   │       ├── bigwigCompare
    │   │   │       │   ├── duplicates_marked
    │   │   │       │   │   └── TSS
    │   │   │       │   └── duplicates_removed
    │   │   │       │       └── TSS
    │   │   │       ├── duplicates_marked
    │   │   │       │   └── TSS
    │   │   │       └── duplicates_removed
    │   │   │           └── TSS
    │   │   ├── multiBamSummary
    │   │   │   ├── duplicates_marked
    │   │   │   └── duplicates_removed
    │   │   ├── plotCorrelation
    │   │   │   ├── duplicates_marked
    │   │   │   └── duplicates_removed
    │   │   ├── plotFingerprint
    │   │   │   ├── duplicates_marked
    │   │   │   └── duplicates_removed
    │   │   ├── plotPCA
    │   │   │   ├── duplicates_marked
    │   │   │   └── duplicates_removed
    │   │   └── plotProfile
    │   │       ├── reference-point
    │   │       │   ├── duplicates_marked
    │   │       │   │   └── TSS
    │   │       │   └── duplicates_removed
    │   │       │       └── TSS
    │   │       └── scale-regions
    │   │           ├── duplicates_marked
    │   │           │   └── TSS
    │   │           └── duplicates_removed
    │   │               └── TSS
    │   ├── macs2
    │   │   ├── H2AZ_vs_Input_TGFb
    │   │   └── H2AZ_vs_Input_WT
    │   │       ├── pooled
    │   │       ├── replicate1
    │   │       └── replicate2
    │   ├── R_Analysis
    │   └── samtools
    │       └── merge
    │           ├── duplicates_marked
    │           └── duplicates_removed
    └── trimmed_data
```
