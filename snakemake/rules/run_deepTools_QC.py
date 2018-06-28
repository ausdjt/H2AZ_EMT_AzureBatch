__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-01-11"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#from snakemake.exceptions import MissingInputException
#import os

"""
Rules for running deepTools QC/QC on ChIP-Seq data
For usage, include this in your workflow.
"""

rule multiBamSummary:
    version:
        0.2
    params:
        dthreads = config["program_parameters"]["deepTools"]["threads"],
        deepTools_dir = config["deepTools_dir"],
        binSize = config["program_parameters"]["deepTools"]["binSize"],
        labels = lambda wildcards: ' '.join("{!s}".format(key) for (key) in config["samples"][wildcards.assayID][wildcards.runID].keys())
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{sample}.Q{qual}.sorted.{suffix}",
               assayID = ASSAY,
               runID = RUNID,
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               sample = config["samples"][ASSAY][RUNID],
               qual = config["alignment_quality"],
               suffix = "bam")
    output:
        npz = "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/multiBamSummary/{duplicates}/results.npz"
    shell:
        """
            {params.deepTools_dir}/multiBamSummary bins --bamfiles {input} \
                                                        --labels {params.labels} \
                                                        --numberOfProcessors {params.dthreads} \
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
                                                   --plotTitle "{params.plotTitle}" \
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
                                           --plotTitle "{params.plotTitle}"
        """

rule bamPEFragmentSize:
    params:
        deepTools_dir = config["deepTools_dir"],
        labels = lambda wildcards: ' '.join("{!s}".format(key) for (key) in config["samples"][wildcards.assayID][wildcards.runID].keys()),
        dthreads = config["program_parameters"]["deepTools"]["threads"]
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{sample}.Q{qual}.sorted.{suffix}",
               assayID = ASSAY,
               runID = RUNID,
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               sample = config["samples"][ASSAY][RUNID],
               qual = config["alignment_quality"],
               suffix = "bam")
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/bamPEFragmentSize/{duplicates}/histogram_duplicates_marked.png"
    shell:
        """
            {params.deepTools_dir}/bamPEFragmentSize --bamfiles {input} \
                                                     --samplesLabel {params.labels} \
                                                     --numberOfProcessors {params.dthreads} \
                                                     --histogram {output}
        """

rule plotFingerprint:
    params:
        deepTools_dir = config["deepTools_dir"],
        plotTitle = lambda wildcards: "BAM PE " + wildcards.duplicates + " fingerprint",
        labels = lambda wildcards: ' '.join("{!s}".format(key) for (key) in config["samples"][wildcards.assayID][wildcards.runID].keys()),
        dthreads = config["program_parameters"]["deepTools"]["threads"]

    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{sample}.Q{qual}.sorted.{suffix}",
               assayID = ASSAY,
               runID = RUNID,
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               sample = config["samples"][ASSAY][RUNID],
               qual = config["alignment_quality"],
               suffix = "bam")
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotFingerprint/{duplicates}/fingerprints_duplicates_marked.png"
    shell:
        """
            {params.deepTools_dir}/plotFingerprint --bamfiles {input} \
                                                   --numberOfProcessors {params.dthreads} \
                                                   --centerReads \
                                                   --plotTitle "{params.plotTitle}" \
                                                   --labels {params.labels} \
                                                   --skipZeros \
                                                   --plotFile {output}
        """
