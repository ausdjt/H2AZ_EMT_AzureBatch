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
    all, run_kallisto, run_STAR, run_htseq, run_cutadapt

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
include:
    include_prefix + "bam_processing_stage2.py"


rule deepTools_QC:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotCorrelation/{duplicates}/heatmap_SpearmanCorr_readCounts.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               suffix = ["png", "tab"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotPCA/{duplicates}/PCA_readCounts.png",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_marked", "duplicates_removed"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotFingerprint/{duplicates}/fingerprints_{duplicates}.png",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_marked", "duplicates_removed"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/bamPEFragmentSize/{duplicates}/histogram_{duplicates}.png",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               duplicates = ["duplicates_marked", "duplicates_removed"])

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{tool}/{duplicates}/{unit}.Q{qual}.sorted.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               tool = "bowtie2",
               duplicates = ["duplicates_marked", "duplicates_removed"],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = ["bam", "bam.bai"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{region}.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               application = "deepTools",
               tool = "plotProfile",
               command = ["reference-point", "scale-regions"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               referencePoint = "TSS",
               plotType = "se",
               mode = ["MNase", "normal"],
               region = "allGenes",
               suffix = ["pdf", "data", "bed"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{command}/{duplicates}/{sample_group}.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               application = "samtools",
               command = "merge",
               duplicates = ["duplicates_marked", "duplicates_removed"],
               sample_group = ["H2AZ-TGFb", "H2AZ-WT", "Input-TGFb", "Input-WT"],
               suffix = ["bam", "bam.bai"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{scaleFactors}/{treatment}_vs_{control}_{mode}_{ratio}_RPKM.bw",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               application = "deepTools",
               tool = "bamCompare",
               mode = ["normal"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               scaleFactors = ["readCount", "SES"],
               treatment = "H2AZ-WT",
               control = "Input-WT"),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{scaleFactors}/{treatment}_vs_{control}_{mode}_{ratio}_RPKM.bw",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               application = "deepTools",
               tool = "bamCompare",
               mode = ["normal"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               scaleFactors = ["readCount", "SES"],
               treatment = "H2AZ-TGFb",
               control = "Input-TGFb")
