__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-29"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for processing SAM/BAM files
For usage, include this in your workflow.
"""

# import other packages
#import os
#import fnmatch
#from snakemake.exceptions import MissingInputException

# set some local variables
#home = os.environ['HOME']

#rule:    version: 0.1

# rules
rule bam_quality_filter:
    params:
        qual = config["alignment_quality"]
    input:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{sample}.bam"
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/quality_filtered/{sample}.Q{qual}.bam"
    shell:
        "/home/apps/samtools/samtools view -b -h -q {params.qual} {input} > {output}"

rule bam_sort:
    params:
        qual = config["alignment_quality"],
        threads = "2"
    input:
        rules.bam_quality_filter.output
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/sorted/{sample}.Q{qual}.sorted.bam"
    shell:
        "/home/apps/samtools/samtools sort -@ {params.threads} {input} -T {wildcards.sample}.Q{params.qual}.sorted -o {output}"

rule bam_mark_duplicates:
    params:
        qual = config["alignment_quality"],
        picard = config["picard"],
        temp = home + config["temp_dir"]
    input:
        rules.bam_sort.output
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{sample}.Q{qual}.sorted.bam"
    shell:
        """
            java -Djava.io.tmpdir={params.temp} \
            -Xmx24G \
            -jar {params.picard} MarkDuplicates \
            INPUT={input}\
            OUTPUT={output}\
            ASSUME_SORTED=TRUE\
            METRICS_FILE={output}.metrics.txt
        """

rule bam_index:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bam_mark_duplicates.output
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{sample}.Q{qual}.sorted.bam.bai"
    shell:
        "/home/apps/samtools/samtools index {input} {output}"

rule bam_rmdup:
    input:
        rules.bam_mark_duplicates.output
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{sample}.Q{qual}.sorted.bam"
    shell:
        "cat {input}> {output}"

rule bam_rmdup_index:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bam_rmdup.output
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{sample}.Q{qual}.sorted.bam.bai"
    shell:
        "/home/apps/samtools/samtools index {input} {output}"
