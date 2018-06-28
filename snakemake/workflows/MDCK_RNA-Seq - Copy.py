_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"

from snakemake.exceptions import MissingInputException
import os
import pdb

#rule:    version: 0.3

#localrules: all, cutadapt_pe, kallisto_quant, star_align_full

home = config["home"]
WORKFLOWDIR = config["WORKFLOWDIR"]
wrapper_dir = home + WORKFLOWDIR + "/snakemake-wrappers/bio"
include_prefix= home + WORKFLOWDIR + "/H2AZ_EMT/snakemake/rules/"
ASSAY = config["ASSAY"] #RNA-Seq
RUNID = config["RUNID"] #NB501086_0082_RDomaschenz_JCSMR_mRNAseq
OUTDIR = config["processed_dir"]

rule all:
    input:
        expand("{assayID}/{runID}/{outDIR}/trimmed_data/{sample}_{suffix}.QT.CA.fastq.gz",
               assayID = ASSAY,
               runID = RUNID,
               sample = config["samples"][ASSAY][RUNID],
               suffix = ["R1_001", "R2_001"],
               outDIR = OUTDIR),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/kallisto/{unit}",
               assayID = ASSAY,
               runID = RUNID,
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"],
               unit = config["samples"][ASSAY][RUNID]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/HTSeq/count/{unit}.txt",
               assayID = ASSAY,
               runID = RUNID,
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"],
               unit = config["samples"][ASSAY][RUNID])

rule cutadapt_pe:
    params:
        trim_params = "-a AGATCGGAAGAGC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --minimum-length=30",
        raw_data = "fastq",
        cutadapt_dir = config["cutadapt_dir"]
    input:
        read1 = lambda wildcards: wildcards.assayID + "/" + wildcards.runID + "/fastq/" + config["samples"][wildcards.assayID][wildcards.runID][wildcards.sample][0],
        read2 = lambda wildcards: wildcards.assayID + "/" + wildcards.runID + "/fastq/" + config["samples"][wildcards.assayID][wildcards.runID][wildcards.sample][1]
    output:
        trimmed_read1 = "{assayID}/{runID}/{outdir}/trimmed_data/{sample}_R1_001.QT.CA.fastq.gz",
        trimmed_read2 = "{assayID}/{runID}/{outdir}/trimmed_data/{sample}_R2_001.QT.CA.fastq.gz"
    run:
        #pdb.set_trace()
        shell("""
            {params.cutadapt_dir}/cutadapt {params.trim_params} \
                                            -o {output.trimmed_read1} \
                                            -p {output.trimmed_read2} \
                                            {input.read1} \
                                            {input.read2}
        """)


rule kallisto_quant:
    params:
        bootstraps = config["program_parameters"]["kallisto"]["bootstraps"],
        threads = 4,
        trim_dir = config["trim_dir"],
        kallisto_dir=config["kallisto_dir"]
    input:
        read1 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R1_001.QT.CA.fastq.gz",
        read2 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R2_001.QT.CA.fastq.gz",
        ki = lambda wildcards: config["references"]["CanFam3.1"]["kallisto"][wildcards.reference_version]
    output:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/kallisto/{unit}"
    shell:
        """
            {params.kallisto_dir}/kallisto quant --index={input.ki} \
                           --output-dir={output} \
                           --threads=4 \
                           --bootstrap-samples={params.bootstraps} \
                           {input.read1} {input.read2}
        """

rule star_align_full:
    version:
        0.4
    params:
        runThreadN = config["program_parameters"]["STAR"]["runThreadN"],
        trim_dir = config["trim_dir"],
        tmp = temp("{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/tmp/{unit}/"),
        star_dir = config["star_dir"]
    input:
        read1 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R1_001.QT.CA.fastq.gz",
        read2 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R2_001.QT.CA.fastq.gz",
        index = lambda wildcards: config["references"]["CanFam3.1"]["STAR"][wildcards.reference_version]
    output:
        bam = "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam"
    shell:
        """
            {params.star_dir}/STAR --runMode alignReads \
                 --runThreadN {params.runThreadN} \
                 --genomeDir {input.index} \
                 --readFilesIn {input.read1} {input.read2}\
                 --readFilesCommand zcat \
                 --outTmpDir {params.tmp} \
				 --outSAMmode Full \
                 --outSAMattributes Standard \
                 --outSAMtype BAM SortedByCoordinate \
                 --outStd BAM_SortedByCoordinate \
                 --alignEndsType EndToEnd\
                 > {output.bam}
        """

rule bam_index_STAR_output:
    version:
        0.2
    input:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam"
    output:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam.bai"
    wrapper:
        "file://" + wrapper_dir + "/samtools/index/wrapper.py"

rule run_htseq_count:
    version:
        0.3
    params:
        htseq_dir = config["HTSeq_dir"]
    input:
        bam = "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam",
        index = "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam.bai",
        gtf = lambda wildcards: config["references"]["CanFam3.1"]["GTF"][wildcards.reference_version]
    output:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/HTSeq/count/{unit}.txt"
    shell:
        """
            {params.htseq_dir}/htseq-count --format=bam \
                                           --order=pos \
                                           --stranded=reverse \
                                           --type=exon \
                                           --idattr=gene_id \
                                           --order=pos \
                                           {input.bam} \
                                           {input.gtf} \
                                           > {output}
        """

# includes for the actual scripts
#include:
#    include_prefix + "perform_cutadapt.py"
#include:
#    include_prefix + "run_kallisto.py"
#include:
#    include_prefix + "run_STAR.py"
