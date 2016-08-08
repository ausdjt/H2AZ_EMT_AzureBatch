__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-06"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for processing SAM/BAM files
For usage, include this in your workflow.
"""

#configfile: ~/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/snakemake/configs/config.json

# import other packages
import os
import fnmatch

# rules
rule all:
 input:
     expand("./H2AZ/processed_data/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.bam", unit = config["Chip"], qual = config["alignment_quality"]),
     expand("./H2AZ/processed_data/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.bam.bai", unit = config["Chip"], qual = config["alignment_quality"]),
     expand("./H2AZ/processed_data/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.bam", unit = config["Chip"], qual = config["alignment_quality"]),
     expand("./H2AZ/processed_data/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.bam.bai", unit = config["Chip"], qual = config["alignment_quality"]),
     expand("./Input/processed_data/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.bam", unit = config["Input"], qual = config["alignment_quality"]),
     expand("./Input/processed_data/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.bam.bai", unit = config["Input"], qual = config["alignment_quality"]),
     expand("./Input/processed_data/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.bam", unit = config["Input"], qual = config["alignment_quality"]),
     expand("./Input/processed_data/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.bam.bai", unit = config["Input"], qual = config["alignment_quality"])
    #  expand("./processed_data/duplicates_removed/{sample}.Q{qual}.sorted.DeDup.bam", sample = ("Input", "H2AZ", "H2A_Bbd"), qual = config["alignment_quality"])

rule bam_rmdup:
    input:
        "./{Type}/processed_data/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.bam"
    output:
        "./{Type}/processed_data/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.bam",
        "./{Type}/processed_data/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.bam.bai"
    shell:
        "samtools rmdup {input} {output[0]}; samtools index {output[0]} {output[1]}"
