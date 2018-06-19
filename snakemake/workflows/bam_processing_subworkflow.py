__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-06-19"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.1

localrules:
    all

home = os.environ['HOME']

wrapper_dir = home + "/Development/snakemake-wrappers/bio"

include_prefix= home + "/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/snakemake/rules/"

configfile: home + "/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/snakemake/configs/config.json"

include:
    include_prefix + "bam_processing.py"

# run parameters as variables
ASSAYID = "ChIP-Seq"
RUNID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"
OUTDIR = "processed_data"
REFVERSION = config["references"]["CanFam3.1"]["version"][0]
QUAL = config["alignment_quality"]

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/bowtie2/{duplicates}/{sample}.Q{qual}.sorted.bam.bai",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               duplicates = ["duplicates_removed", "duplicates_marked"]
               sample = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"])
