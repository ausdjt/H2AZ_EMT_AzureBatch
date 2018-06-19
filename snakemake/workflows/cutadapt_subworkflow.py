__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-06-18"

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
    include_prefix + "perform_cutadapt.py"

# run parameters as variables
ASSAYID = "ChIP-Seq"
RUNID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"
OUTDIR = "processed_data"

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{trim_data}/{unit}_{suffix}.QT.CA.fastq.gz",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               trim_data = config["trim_dir"],
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
               suffix = ["R1_001", "R2_001"])
