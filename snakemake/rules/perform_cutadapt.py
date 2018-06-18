__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-27"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

# this set of rules is meant to be imported by the master workflow document

from snakemake.exceptions import MissingInputException
import os
import pdb
home = os.environ['HOME']

"""
Rules for trimming reads with cutadapt
(http://cutadapt.readthedocs.org/en/latest/guide.html#illumina-truseq)

For usage, include this in your workflow.
"""

rule cutadapt_pe:
    params:
        trim_params = "-a AGATCGGAAGAGC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --minimum-length=30",
        raw_data = "fastq",
        cutadapt_dir = home + "miniconda3/envs/snakemake/bin"
    input:
        read1 = lambda wildcards: wildcards.assayID + "/" + wildcards.runID + "/fastq/" + config["samples"][wildcards.assayID][wildcards.runID][wildcards.sample][0],
        read2 = lambda wildcards: wildcards.assayID + "/" + wildcards.runID + "/fastq/" + config["samples"][wildcards.assayID][wildcards.runID][wildcards.sample][1]
    output:
        trimmed_read1 = "{assayID}/{runID}/{outdir}/trimmed_data/{sample}_R1_001.QT.CA.fastq.gz",
        trimmed_read2 = "{assayID}/{runID}/{outdir}/trimmed_data/{sample}_R2_001.QT.CA.fastq.gz"
    run:
        pdb.set_trace()
        shell("""
            {params.cutadapt_dir}/cutadapt {params.trim_params} \
                                            -o {output.trimmed_read1} \
                                            -p {output.trimmed_read2} \
                                            {input.read1} \
                                            {input.read2}
        """)
