__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-06-19"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os

"""
Rules for running deepTools QC/QC on ChIP-Seq data
For usage, include this in your workflow.
"""
rule:
    version: 0.1

localrules:
    all

# run parameters as variables
ASSAY = config["ASSAY"]
RUNID = config["RUNID"]
OUTDIR = config["processed_dir"]
REFVERSION = config["references"]["CanFam3.1"]["version"][0]
QUAL = config["alignment_quality"]
home = os.environ['HOME']
WORKFLOWDIR = config["WORKFLOWDIR"]
wrapper_dir = home + WORKFLOWDIR + "snakemake-wrappers/bio"
include_prefix= home + WORKFLOWDIR + "H2AZ_EMT/snakemake/rules/"
configfile: home + WORKFLOWDIR + "H2AZ_EMT/snakemake/configs/config.json"

# includes for the actual scripts
include:
    include_prefix + "run_deepTools_QC.py"

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotFingerprint/{duplicates}/fingerprints_duplicates_marked.png",
               assayID = ASSAY,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               duplicates = ["duplicates_removed", "duplicates_marked"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotCorrelation/{duplicates}/heatmap_SpearmanCorr_readCounts.{suffix}",
               assayID = ASSAY,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               duplicates = ["duplicates_removed", "duplicates_marked"],
               suffix = ["png", "tab"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotPCA/{duplicates}/PCA_readCounts.png",
               assayID = ASSAY,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               duplicates = ["duplicates_removed", "duplicates_marked"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/bamPEFragmentSize/{duplicates}/histogram_duplicates_marked.png",
               assayID = ASSAY,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               duplicates = ["duplicates_removed", "duplicates_marked"])
