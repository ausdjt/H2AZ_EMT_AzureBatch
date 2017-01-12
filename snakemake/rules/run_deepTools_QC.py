__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-01-11"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os

"""
Rules for running deepTools QC/QC on ChIP-Seq data
For usage, include this in your workflow.
"""

rule multiBamSummary:
    version:
        0.2
    params:
        deepTools_dir = config["deepTools_dir"],
        binSize = config["program_parameters"]["deepTools"]["binSize"],
        labels = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"]
    threads:
        24
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{unit}.Q{qual}.sorted.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_marked"],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = "MkDup.bam")
    output:
        npz = "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/multiBamSummary/duplicates_marked/results.npz"
    shell:
        """
            {params.deepTools_dir}/multiBamSummary bins --bamfiles {input} \
                                                        --labels {params.labels} \
                                                        --numberOfProcessors {threads} \
                                                        --centerReads \
                                                        --binSize {params.binSize} \
                                                        --outFileName {output.npz}
        """

rule multiBamSummary_deduplicated:
    version:
        0.2
    params:
        deepTools_dir = config["deepTools_dir"],
        binSize = config["program_parameters"]["deepTools"]["binSize"],
        labels = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"]
    threads:
        24
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{unit}.Q{qual}.sorted.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_removed"],
               qual = config["alignment_quality"],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               suffix = "DeDup.bam")
    output:
        npz = "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/multiBamSummary/duplicates_removed/results.npz"
    shell:
        """
            {params.deepTools_dir}/multiBamSummary bins --bamfiles {input} \
                                                        --labels {params.labels} \
                                                        --numberOfProcessors {threads} \
                                                        --centerReads \
                                                        --binSize {params.binSize} \
                                                        --outFileName {output.npz}
        """


rule plotCorrelation_heatmap:
    params:
        deepTools_dir = config["deepTools_dir"],
        plotTitle = lambda wildcards: "Correlation heatmap - " + wildcards.duplicates
    input:
        npz = "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/multiBamSummary/{duplicates}/results.npz"
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotCorrelation/{duplicates}/heatmap_SpearmanCorr_readCounts.png",
        "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotCorrelation/{duplicates}/heatmap_SpearmanCorr_readCounts.tab"
    shell:
        """
            {params.deepTools_dir}/plotCorrelation -in {input.npz} \
                                                   --corMethod spearman \
                                                   --skipZeros \
                                                   --plotTitle {params.plotTitle} \
                                                   --whatToPlot heatmap \
                                                   --colorMap RdYlBu \
                                                   --plotNumbers \
                                                   -o {output[0]} \
                                                   --outFileCorMatrix {output[1]}
        """

rule plotPCA:
    params:
        deepTools_dir = config["deepTools_dir"],
        plotTitle = lambda wildcards: "PCA - " + wildcards.duplicates
    input:
        npz = "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/multiBamSummary/{duplicates}/results.npz"
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotPCA/{duplicates}/PCA_readCounts.png"
    shell:
        """
            {params.deepTools_dir}/plotPCA -in {input.npz} \
                                           -o {output} \
                                           -plotTitle {params.plotTitle}
        """

rule bamPEFragmentSize:
    params:
        deepTools_dir = config["deepTools_dir"],
        labels = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"]
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{unit}.Q{qual}.sorted.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_marked"],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = "MkDup.bam")
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/bamPEFragmentSize/duplicates_marked/histogram.png"
    shell:
        """
            {params.deepTools_dir}/bamPEFragmentSize --bamfiles {input} \
                                                     --samplesLabel {params.labels} \
                                                     --numberOfProcessors {threads} \
                                                     --histogram {output}
        """

rule bamPEFragmentSize_deduplicated:
    params:
        deepTools_dir = config["deepTools_dir"],
        labels = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
        plotTitle = lambda wildcards: "BAM PE " + wildcards.duplicates + " fragment size"
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{unit}.Q{qual}.sorted.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_removed"],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = "DeDup.bam")
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/bamPEFragmentSize/duplicates_removed/histogram.png"
    shell:
        """
            {params.deepTools_dir}/bamPEFragmentSize --bamfiles {input} \
                                                     --samplesLabel {params.labels} \
                                                     --numberOfProcessors {threads} \
                                                     --plotTitle {params.plotTitle} \
                                                     --histogram {output}
        """

rule plotFingerprint:
    params:
        deepTools_dir = config["deepTools_dir"],
        plotTitle = lambda wildcards: "BAM PE " + wildcards.duplicates + " fingerprint"
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_marked"],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = "MkDup.bam")
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotFingerprint/{duplicates}/fingerprints.png"
    shell:
        """
            {params.deepTools_dir}/plotFingerprint --bamfiles {input} \
                                                   --numberOfProcessors {threads} \
                                                   --centerReads \
                                                   --plotTitle {params.plotTitle} \
                                                   --skipZeros \
                                                   --plotFile {output}
        """

rule plotFingerprint_deduplicated:
    params:
        deepTools_dir = config["deepTools_dir"],
        plotTitle = lambda wildcards: "BAM PE " + wildcards.duplicates + " fingerprint"
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_removed"],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = "DeDup.bam")
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotFingerprint/{duplicates}/fingerprints.png"
    shell:
        """
            {params.deepTools_dir}/plotFingerprint --bamfiles {input} \
                                                   --numberOfProcessors {threads} \
                                                   --centerReads \
                                                   --plotTitle {params.plotTitle} \
                                                   --skipZeros \
                                                   --plotFile {output}
        """
