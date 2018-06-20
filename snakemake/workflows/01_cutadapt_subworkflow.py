__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-06-18"

from snakemake.exceptions import MissingInputException
import os

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
    include_prefix + "perform_cutadapt.py"

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{trim_data}/{unit}_{suffix}.QT.CA.fastq.gz",
               assayID = ASSAY,
               runID = RUNID,
               outdir = OUTDIR,
               trim_data = config["trim_dir"],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               suffix = ["R1_001", "R2_001"])
