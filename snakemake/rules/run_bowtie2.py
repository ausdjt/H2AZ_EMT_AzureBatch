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

rule bowtie2_pe:
    params:
        threads= config["program_parameters"]["bt2_params"]["threads"],
        max_in= config["program_parameters"]["bt2_params"]["max_insert"],
        bt2_index= home + config["references"]["CanFam3.1"]["genome"]
    input:
        read1="{assayID}/{runID}/{outdir}/trimmed_data/{unit}_R1_001.QT.CA.fastq.gz",
        read2="{assayID}/{runID}/{outdir}/trimmed_data/{unit}_R2_001.QT.CA.fastq.gz"
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{unit}.bam")
    shell:
        """
            bowtie2 \
            -x {params.bt2_index}\
            --no-mixed \
            --no-discordant \
            --maxins {params.max_in} \
            --threads {params.threads}\
            --rg-id '{wildcards.unit}' \
            --rg 'LB:{wildcards.unit}' \
            --rg 'SM:{wildcards.unit}' \
            --rg 'PL:Illumina' \
            --rg 'PU:NA' \
            -1 {input.read1} \
            -2 {input.read2} \
            | samtools view -Sb - > {output}
        """

rule bam_stats:
    input:
        rules.bowtie2_pe.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{unit}.bam.stats.txt")
    shell:
        """
            samtools flagstat {input} > {output}
        """

rule bam_extract_unmapped_reads:
    input:
        rules.bowtie2_pe.output
    output:
        temp("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/unmapped_reads/{unit}_unmapped.bam")
    shell:
        "samtools view -f 4 -b {input} > {output}"

rule bam_sort_unmapped_reads:
    input:
        rules.bam_extract_unmapped_reads.output
    output:
        temp("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/unmapped_reads/sorted/{unit}_unmapped.sorted.bam")
    shell:
        "samtools sort -n {input} -T {wildcards.unit}.sorted -o {output}"

rule unmapped_reads_to_pe_fastq:
    version:
        "0.1"
    input:
        rules.bam_sort_unmapped_reads.output
    output:
        temp("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/unmapped_reads/fastq/{unit}_unmapped_r1.fastq"),
        temp("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/unmapped_reads/fastq/{unit}_unmapped_r2.fastq")
    shell:
        """
            bedtools bamtofastq -i {input} \
                                -fq {output[0]} \
                                -fq2 {output[1]}
        """

rule gzip_unmapped_fastq:
    version:
        "0.1"
    input:
        rules.unmapped_reads_to_pe_fastq.output
    output:
        temp("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/unmapped_reads/fastq/{unit}_unmapped_r1.fastq.gz"),
        temp("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/unmapped_reads/fastq/{unit}_unmapped_r2.fastq.gz")
    shell:
        "gzip {input[0]}; gzip {input[1]}"

rule cutadapt_pe_unmapped:
    """Trims given paired-end reads with given parameters"""
    version:
        "0.2"
    params:
        trim_params = config["program_parameters"]["cutadapt"]["trim_params"],
        cutadapt_dir = home + config["cutadapt_dir"]
    input:
        rules.gzip_unmapped_fastq.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/unmapped_reads/fastq/{unit}_r1.QT.CA.fastq.gz"),
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/unmapped_reads/fastq/{unit}_r2.QT.CA.fastq.gz")
    shell:
        """
            {params.cutadapt_dir}/cutadapt {params.trim_params} \
            -o {output[0]} -p {output[1]} \
            {input[0]} {input[1]}
        """

rule bowtie2_pe_unmapped_reads:
    version:
        "0.1"
    message:
        "Running second round of bowtie2 alignments on previously unmapped reads..."
    params:
        threads = config["program_parameters"]["bt2_params"]["threads"],
        max_in = config["program_parameters"]["bt2_params"]["max_insert"],
        bt2_index = home + config["references"]["CanFam3.1"]["genome"]
    input:
        rules.cutadapt_pe_unmapped.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/second_alignment/{unit}.2ndAlignment.bam")
    shell:
        """
            bowtie2 \
            -x {params.bt2_index}\
            --maxins {params.max_in} \
            --threads {params.threads}\
            --very-sensitive\
            --rg-id '{wildcards.unit}' \
            --rg 'LB:{wildcards.unit}' \
            --rg 'SM:{wildcards.unit}' \
            --rg 'PL:Illumina' \
            --rg 'PU:NA' \
            -1 {input[0]} \
            -2 {input[1]} \
            | samtools view -Sb - > {output}
        """

rule bam_stats_2nd_alignment:
    input:
        rules.bowtie2_pe_unmapped_reads.output
    output:
        protected("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/second_alignment/{unit}.2ndAlignment.bam.stats.txt")
    shell:
        """
            samtools flagstat {input} > {output}
        """


# rule bowtie2_pe_multi_mapper:
#     params:
#         threads = config["bt2_params"]["threads"],
#         max_in = config["bt2_params"]["max_insert"],
#         bt2_index = config["references"]["genome"]
#     input:
#         "NB501086_0034_MNekrasov_JCSMR_ChIPseq/{unit}_R1_001.fastq.gz",
#         "NB501086_0034_MNekrasov_JCSMR_ChIPseq/{unit}_R2_001.fastq.gz"
#     output:
#         "NB501086_0034_MNekrasov_JCSMR_ChIPseq/processed_data/multi_mapped/{unit}.bam"
#     shell:
#         """
#             bowtie2 \
#             -x {params.bt2_index}\
#             --no-mixed \
#             --no-discordant \
#             --maxins {params.max_in} \
#             --threads {params.threads}\
#             -k 60\
#             --rg-id '{wildcards.unit}' \
#             --rg 'LB:{wildcards.unit}' \
#             --rg 'SM:{wildcards.unit}' \
#             --rg 'PL:Illumina' \
#             --rg 'PU:NA' \
#             -1 {input[0]} \
#             -2 {input[1]} \
#             | samtools view -Sb - > {output}
#         """
