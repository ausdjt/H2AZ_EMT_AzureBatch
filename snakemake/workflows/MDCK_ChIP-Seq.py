_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-12-05"

from snakemake.exceptions import MissingInputException
import snakemake.utils
import os
import pdb

rule:
    version: 0.1

localrules:
    all

home = os.environ['HOME']

wrapper_dir = home + "/Development/snakemake-wrappers/bio"

include_prefix = home + "/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/snakemake/rules/"

#include:
#    include_prefix + "perform_cutadapt.py"
include:
    include_prefix + "run_bowtie2.py"
include:
    include_prefix + "bam_processing.py"
include:
    include_prefix + "run_deepTools_QC.py"

include:
    include_prefix + "run_deepTools.py"

# run parameters as variables
RUNID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"
ASSAYID = "ChIP-Seq"
OUTDIR = config["processed_dir"]
REFVERSION = config["references"]["CanFam3.1"]["version"][0]
QUALITY = config["alignment_quality"]

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{tool}/{duplicates}/{unit}.Q{qual}.sorted.{suffix}",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               tool = "bowtie2",
               duplicates = ["duplicates_marked", "duplicates_removed"],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               qual = QUALITY,
               suffix = ["bam", "bam.bai"])

rule deepTools_QC:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotCorrelation/{duplicates}/heatmap_SpearmanCorr_readCounts.{suffix}",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               duplicates = ["duplicates_marked", "duplicates_removed"],
               suffix = ["png", "tab"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotPCA/{duplicates}/PCA_readCounts.png",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               duplicates = ["duplicates_marked", "duplicates_removed"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotFingerprint/{duplicates}/fingerprints_{duplicates}.png",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               duplicates = ["duplicates_marked", "duplicates_removed"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/bamPEFragmentSize/{duplicates}/histogram_{duplicates}.png",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               duplicates = ["duplicates_marked", "duplicates_removed"])

rule run_computeMatrix:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{region}_{mode}.matrix.gz",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               application = "deepTools",
               tool = "computeMatrix",
               command = ["reference-point", "scale-regions"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               referencePoint = "TSS",
               region = ["allGenes", "Tan_EMT_up", "Tan_EMT_down"],
               mode = ["MNase", "normal"])

rule allGenes_plots:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{region}.{suffix}",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               application = "deepTools",
               tool = "plotProfile",
               command = ["reference-point", "scale-regions"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               referencePoint = "TSS",
               plotType = "se",
               mode = ["MNase", "normal"],
               region = "allGenes",
               suffix = ["pdf", "data", "bed"])

rule bamCoverage_replicates:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/merged_replicates/{sampleGroup}_{mode}_{norm}.bw",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               application = "deepTools",
               tool = "bamCoverage",
               mode = ["normal", "MNase"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               sampleGroup = ["H2AZ-TGFb", "H2AZ-WT", "Input-TGFb", "Input-WT"],
               norm = "RPKM")

rule bigwigCompare_replicates:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{scaleFactors}/{treatment}_vs_{control}_{mode}_{ratio}_{norm}.bw",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               application = "deepTools",
               tool = "bigwigCompare",
               mode = ["normal"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               scaleFactors = ["readCount", "SES"],
               treatment = "H2AZ-WT",
               control = "Input-WT",
               ratio = ["log2", "subtract"],
               norm = "RPKM"),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{scaleFactors}/{treatment}_vs_{control}_{mode}_{ratio}_{norm}.bw",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               application = "deepTools",
               tool = "bigwigCompare",
               mode = ["normal"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               scaleFactors = ["readCount", "SES"],
               treatment = "H2AZ-TGFb",
               control = "Input-TGFb",
               ratio = ["log2", "subtract"],
               norm = "RPKM"),
