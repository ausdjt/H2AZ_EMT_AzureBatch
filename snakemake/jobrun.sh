#!/usr/bin/env bash
cd /data
#snakemake --latency-wait 60 --snakefile ./Development/H2AZ_EMT/snakemake/workflows/MDCK_ChIP-Seq.py --configfile ./Development/H2AZ_EMT/snakemake/configs/config.json --config ASSAY=ChIP-Seq RUNID=NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq WORKFLOWDIR=Development --jobs -pr
snakemake --latency-wait 60 --snakefile ./Development/H2AZ_EMT/snakemake/workflows/MDCK_RNA-Seq.py --configfile ./Development/H2AZ_EMT/snakemake/configs/config.json --config ASSAY=RNA-Seq RUNID=NB501086_0082_RDomaschenz_JCSMR_mRNAseq WORKFLOWDIR=Development --jobs -pr
