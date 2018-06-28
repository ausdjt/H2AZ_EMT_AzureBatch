__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-12-05"

from snakemake.exceptions import MissingInputException
import snakemake.utils
import os
import pdb

home = config["home"]
WORKFLOWDIR = config["WORKFLOWDIR"]
wrapper_dir = home + WORKFLOWDIR + "/snakemake-wrappers/bio"
include_prefix= home + WORKFLOWDIR + "/H2AZ_EMT/snakemake/rules/"
ASSAY = config["ASSAY"] #ChIP-Seq
RUNID = config["RUNID"] #NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq
OUTDIR = config["processed_dir"]
REFVERSION = config["references"]["CanFam3.1"]["version"][0]
QUALITY = config["alignment_quality"]

rule all:
    input:
        expand("{assayID}/{runID}/{outDIR}/trimmed_data/{sample}_{suffix}.QT.CA.fastq.gz",
               assayID = ASSAY,
               runID = RUNID,
               sample = config["samples"][ASSAY][RUNID],
               suffix = ["R1_001", "R2_001"],
               outDIR = OUTDIR),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{sample}.Q{qual}.sorted.bam.bai",
               assayID = ASSAY,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               duplicates = ["duplicates_removed", "duplicates_marked"],
               qual=QUALITY,
               sample = config["samples"][ASSAY][RUNID]),
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

include: include_prefix + "perform_cutadapt.py"
include: include_prefix + "bam_processing.py"
include: include_prefix + "run_bowtie2.py"
#include: include_prefix + "run_deepTools.py"
include: include_prefix + "run_deepTools_QC.py"

#include:
#    include_prefix + "pooled_replicates_processing.py"




# rule run_computeMatrix_pooled_replicates:
#     input:
#         expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/{sampleGroup}_{region}_{mode}.matrix.gz",
#                assayID = ASSAYID,
#                runID = RUNID,
#                outdir = OUTDIR,
#                reference_version = REFVERSION,
#                application = "deepTools",
#                tool = "computeMatrix",
#                command = ["reference-point", "scale-regions"],
#                duplicates = ["duplicates_marked", "duplicates_removed"],
#                referencePoint = "TSS",
#                sampleGroup = ["H2AZ-TGFb", "H2AZ-WT", "Input-TGFb", "Input-WT"],
#                region = ["allGenes", "Tan_EMT_up", "Tan_EMT_down"],
#                mode = ["MNase", "normal"])

# rule allGenes_plots:
#     input:
#         expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{region}.{suffix}",
#                assayID = ASSAYID,
#                runID = RUNID,
#                outdir = OUTDIR,
#                reference_version = REFVERSION,
#                application = "deepTools",
#                tool = "plotProfile",
#                command = ["reference-point", "scale-regions"],
#                duplicates = ["duplicates_marked", "duplicates_removed"],
#                referencePoint = "TSS",
#                plotType = "se",
#                mode = ["MNase", "normal"],
#                region = "allGenes",
#                suffix = ["pdf", "data", "bed"])
#
# rule bamCoverage_replicates:
#     input:
#         expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/merged_replicates/{sampleGroup}_{mode}_{norm}.bw",
#                assayID = ASSAYID,
#                runID = RUNID,
#                outdir = OUTDIR,
#                reference_version = REFVERSION,
#                application = "deepTools",
#                tool = "bamCoverage",
#                mode = ["normal", "MNase"],
#                duplicates = ["duplicates_marked", "duplicates_removed"],
#                sampleGroup = ["H2AZ-TGFb", "H2AZ-WT", "Input-TGFb", "Input-WT"],
#                norm = "RPKM")
#
# rule bigwigCompare_replicates:
#     input:
#         expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{scaleFactors}/{treatment}_vs_{control}_{mode}_{ratio}_{norm}.bw",
#                assayID = ASSAYID,
#                runID = RUNID,
#                outdir = OUTDIR,
#                reference_version = REFVERSION,
#                application = "deepTools",
#                tool = "bigwigCompare",
#                mode = ["normal"],
#                duplicates = ["duplicates_marked", "duplicates_removed"],
#                scaleFactors = ["readCount", "SES"],
#                treatment = "H2AZ-WT",
#                control = "Input-WT",
#                ratio = ["log2", "subtract"],
#                norm = "RPKM"),
#         expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{scaleFactors}/{treatment}_vs_{control}_{mode}_{ratio}_{norm}.bw",
#                assayID = ASSAYID,
#                runID = RUNID,
#                outdir = OUTDIR,
#                reference_version = REFVERSION,
#                application = "deepTools",
#                tool = "bigwigCompare",
#                mode = ["normal"],
#                duplicates = ["duplicates_marked", "duplicates_removed"],
#                scaleFactors = ["readCount", "SES"],
#                treatment = "H2AZ-TGFb",
#                control = "Input-TGFb",
#                ratio = ["log2", "subtract"],
#                norm = "RPKM")
