# H2AZ_EMT
Bash and R scripts associated with analysis of H2A.Z ChIP-Seq data from TGFb-induced EMT in MDCK cells. Pre-release.

# some requirements and status of different pipeline stages
## tools & binaries
snakemake [>4.0]
bowtie2 [>2.1]
cutadapt
samtools
deepTools
R/Bioconductor [>3.4]

## pipeline stages
cutadapt_subworkflow.py - OK
bowtie2_subworkflow.py - OK
bam_processing_subworkflow.py - OK
merge_replicates_subworkflow.py - seems OK?
deepTools_QC_subworkflow.py - OK
