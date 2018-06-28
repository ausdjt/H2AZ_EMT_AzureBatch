__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-06-20"

from snakemake.exceptions import MissingInputException
import snakemake.utils
import os

rule:
    version: 0.1

# run parameters as variables
ASSAY = config["ASSAY"]
RUNID = config["RUNID"]
OUTDIR = config["processed_dir"]
REFVERSION = config["references"]["CanFam3.1"]["version"][0]
QUAL = config["alignment_quality"]
home = os.environ['HOME']
WORKFLOWDIR = config["WORKFLOWDIR"]
wrapper_dir = home + WORKFLOWDIR + "snakemake-wrappers/bio"
include_prefix = home + WORKFLOWDIR + "H2AZ_EMT/snakemake/rules/"
configfile: home + WORKFLOWDIR + "H2AZ_EMT/snakemake/configs/config.json"

# includes for the actual scripts
include:
    include_prefix + "run_deepTools_pooled_data.py"

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/{sampleGroup}_{region}_{mode}.matrix.gz",
               assayID = ASSAY,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               application = "deepTools",
               tool = "computeMatrix",
               command = ["reference-point", "scale-regions"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               referencePoint = "TSS",
               sampleGroup = ["H2AZ-TGFb", "H2AZ-WT", "Input-TGFb", "Input-WT"],
               region = ["allGenes", "TanEMTdown", "TanEMTup", "qPCRGenesUp", "qPCRGenesDown", "random100up", "random100down"],
               mode = ["MNase", "normal"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{scaleFactors}/{treatment}_vs_{control}_{mode}_{ratio}_{norm}.bw",
               assayID = ASSAY,
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
               assayID = ASSAY,
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
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/allSamples_{plotType}.{mode}.{region}.{suffix}",
                assayID = ASSAY,
                runID = RUNID,
                outdir = OUTDIR,
                reference_version = REFVERSION,
                application = "deepTools",
                tool = "plotProfile",
                command = ["reference-point", "scale-regions"],
                duplicates = ["duplicates_marked", "duplicates_removed"],
                referencePoint = "TSS",
                plotType = "se",
                region = ["allGenes", "TanEMTup", "TanEMTdown", "qPCRGenesUp", "qPCRGenesDown", "random100up", "random100down"],
                mode = ["MNase", "normal"],
                suffix = ["pdf", "bed", "data"])
