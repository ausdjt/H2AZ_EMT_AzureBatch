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

# set some local variables
home = os.environ['HOME']

def bam_merge_input(wildcards):
    fn = []
    path = "/".join((wildcards["assayID"],
                     wildcards["runID"],
                     wildcards["outdir"],
                     wildcards["reference_version"],
                     "bowtie2",
                     wildcards["duplicates"]))
    for i in config["samples"]["ChIP-Seq"]["replicates"][wildcards["sample_group"]]:
        fn.append("/".join((path, ".".join((i, "Q20.sorted.bam")))))
    return(fn)

rule bam_merge:
    version:
        0.2
    params:
        cwd = os.getcwd()
    threads:
        lambda wildcards: int(str(config["program_parameters"]["bt2_params"]["threads"]).strip("['']"))
    input:
        bam_merge_input
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{command}/{duplicates}/{sample_group}.bam")
    run:
        if (len(input) > 1):
            shell("samtools merge --threads {threads} {output} {input}")
        else:
            shell("ln -s {params.cwd}/{input} {output}")

rule index_merged_bam:
    input:
        rules.bam_merge.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{command}/{duplicates}/{sample_group}.bam.bai")
    shell:
        "samtools index {input} {output}"
