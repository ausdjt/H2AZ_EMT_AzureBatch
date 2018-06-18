__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-26"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for aligning paired-end reads using bowtie2.

For use, include in your workflow.
"""

import os
import fnmatch
from snakemake.exceptions import MissingInputException

# set local variables
home = os.environ['HOME']

subworkflow cutadapt:
    workdir: "."
    snakefile: "cutadapt_subworkflow.smr"

rule bowtie2_pe:
    params:
        threads= config["program_parameters"]["bt2_params"]["threads"],
        max_in= config["program_parameters"]["bt2_params"]["max_insert"],
        bt2_index= home + config["references"]["CanFam3.1"]["genome"]
    input:
        read1=cutadapt("{assayID}/{runID}/{outdir}/trimmed_data/{sample}_R1_001.QT.CA.fastq.gz"),
        read2=cutadapt("{assayID}/{runID}/{outdir}/trimmed_data/{sample}_R2_001.QT.CA.fastq.gz")
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{sample}.bam")
    shell:
        """
            bowtie2 \
            -x {params.bt2_index}\
            --no-mixed \
            --no-discordant \
            --maxins {params.max_in} \
            --threads {params.threads}\
            --rg-id '{wildcards.sample}' \
            --rg 'LB:{wildcards.sample}' \
            --rg 'SM:{wildcards.sample}' \
            --rg 'PL:Illumina' \
            --rg 'PU:NA' \
            -1 {input.read1} \
            -2 {input.read2} \
            | samtools view -Sb - > {output}
        """
