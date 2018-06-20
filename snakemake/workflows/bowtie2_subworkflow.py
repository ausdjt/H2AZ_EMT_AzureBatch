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

include:
    include_prefix + "run_bowtie2.py"

# run parameters as variables
RUNID = config["RUNID"] #"NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"
ASSAYID = config["ASSAYID"] #"ChIP-Seq"
OUTDIR = config["processed_dir"]
REFVERSION = config["references"]["CanFam3.1"]["version"][0]

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{tool}/{unit}.bam",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               tool = "bowtie2",
               unit = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"])
