_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-12-05"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.1

localrules:
    all, run_kallisto, run_STAR, run_htseq, run_cutadapt

home = os.environ['HOME']

wrapper_dir = home + "/Development/snakemake-wrappers/bio"

include_prefix = home + "/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/snakemake/rules/"

# include:
# include:
#     include_prefix + "perform_cutadapt.py"
include:
    include_prefix + "run_bowtie2.py"
include:
    include_prefix + "bam_processing.py"
include:
    include_prefix + "run_deepTools_QC.py"

rule run_cutadapt:
    input:
        expand("{assayID}/{runID}/{outdir}/{trim_data}/{unit}_{suffix}.QT.CA.fastq.gz",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               trim_data = config["trim_dir"],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               suffix = ["R1_001", "R2_001"])

rule deepTools_QC:
    input:
            expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotCorrelation/{duplicates}/heatmap_SpearmanCorr_readCounts.png",
                   assayID = "ChIP-Seq",
                   runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
                   outdir = config["processed_dir"],
                   reference_version = config["references"]["CanFam3.1"]["version"][0],
                   duplicates = ["duplicates_removed", "duplicates_removed"],
                   suffix = ["png", "tab"]),
            expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotPCA/{duplicates}/PCA_readCounts.png",
                   assayID = "ChIP-Seq",
                   runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
                   outdir = config["processed_dir"],
                   reference_version = config["references"]["CanFam3.1"]["version"][0],
                   duplicates = ["duplicates_removed", "duplicates_removed"]),
            expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/bamPEFragmentSize/{duplicates}/histogram.png",
                   assayID = "ChIP-Seq",
                   runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
                   outdir = config["processed_dir"],
                   reference_version = config["references"]["CanFam3.1"]["version"][0],
                   duplicates = "duplicates_marked"),
            expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotFingerprint/{duplicates}/fingerprints.png",
                   assayID = "ChIP-Seq",
                   runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
                   outdir = config["processed_dir"],
                   reference_version = config["references"]["CanFam3.1"]["version"][0],
                   duplicates = "duplicates_marked")

rule deepTools_QC_deduplicated:
    input:
            expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/plotFingerprint/{duplicates}/fingerprints.png",
                   assayID = "ChIP-Seq",
                   runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
                   outdir = config["processed_dir"],
                   reference_version = config["references"]["CanFam3.1"]["version"][0],
                   duplicates = "duplicates_removed")
            expand("{assayID}/{runID}/{outdir}/{reference_version}/deepTools/bamPEFragmentSize/{duplicates}/histogram.png",
                   assayID = "ChIP-Seq",
                   runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
                   outdir = config["processed_dir"],
                   reference_version = config["references"]["CanFam3.1"]["version"][0],
                   duplicates = "duplicates_removed"),


rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_marked/{unit}.Q{qual}.sorted.MkDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = ["bam", "bam.bai"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/duplicates_removed/{unit}.Q{qual}.sorted.DeDup.{suffix}",
               assayID = "ChIP-Seq",
               runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"][0],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               qual = config["alignment_quality"],
               suffix = ["bam", "bam.bai"])
