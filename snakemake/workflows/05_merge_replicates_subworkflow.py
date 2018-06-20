__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-02-01"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for processing SAM/BAM files
For usage, include this in your workflow.
"""

# import other packages
import os
import fnmatch
from snakemake.exceptions import MissingInputException
# run parameters as variables
ASSAY = config["ASSAY"]
RUNID = config["RUNID"]
OUTDIR = config["processed_dir"]
REFVERSION = config["references"]["CanFam3.1"]["version"][0]
QUAL = config["alignment_quality"]
home = os.environ['HOME']
WORKFLOWDIR = config["WORKFLOWDIR"]
wrapper_dir = home + WORKFLOWDIR + "snakemake-wrappers/bio"
include_prefix = home + WORKFLOWDIR + "H2AZ_EMT/snakemake/rules/"

configfile:
    home + WORKFLOWDIR + "H2AZ_EMT/snakemake/configs/config.json"

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{command}/{duplicates}/{sampleGroup}.{suffix}",
               assayID = ASSAY,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               application = "samtools",
               command = "merge",
               duplicates = ["duplicates_marked", "duplicates_removed"],
               sampleGroup = ["H2AZ-TGFb", "H2AZ-WT", "Input-TGFb", "Input-WT"],
               suffix = ["bam.bai"]),


rule bam_merge:
    version:
        0.2
    params:
        cwd = os.getcwd()
    threads:
        lambda wildcards: int(str(config["program_parameters"]["bt2_params"]["threads"]).strip("['']"))
    input:
        lambda wildcards: expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{unit}.Q{qual}.{suffix}",
                                 assayID = wildcards["assayID"],
                                 runID = wildcards["runID"],
                                 outdir = wildcards["outdir"],
                                 reference_version = wildcards["reference_version"],
                                 duplicates = wildcards["duplicates"],
                                 unit = config["samples"]["ChIP-Seq"]["replicates"][wildcards["sampleGroup"]],
                                 qual = config["alignment_quality"],
                                 suffix = "sorted.bam")
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{command}/{duplicates}/{sampleGroup}.bam")
    run:
        if (len(input) > 1):
            shell("samtools merge --threads {threads} {output} {input}")
        else:
            shell("ln -s {params.cwd}/{input} {output}")

rule index_merged_bam:
    input:
        rules.bam_merge.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{command}/{duplicates}/{sampleGroup}.bam.bai")
    shell:
        "samtools index {input} {output}"
