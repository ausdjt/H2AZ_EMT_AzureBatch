_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.3

localrules:
    all, run_kallisto, run_STAR, run_htseq, run_cutadapt

home = os.environ['HOME']

wrapper_dir = home + "/Development/snakemake-wrappers/bio"

include_prefix= home + "/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/snakemake/rules/"

include:
    include_prefix + "perform_cutadapt.py"
include:
    include_prefix + "run_kallisto.py"
include:
    include_prefix + "run_STAR.py"

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/kallisto/{unit}",
               assayID = "RNA-Seq",
               runID = ["NB501086_0082_RDomaschenz_JCSMR_mRNAseq"],
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"],
               unit = config["samples"]["RNA-Seq"]["NB501086_0082_RDomaschenz_JCSMR_mRNAseq"]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/HTSeq/count/{unit}.txt",
               assayID = "RNA-Seq",
               runID = "NB501086_0082_RDomaschenz_JCSMR_mRNAseq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"],
               unit = config["samples"]["RNA-Seq"]["NB501086_0082_RDomaschenz_JCSMR_mRNAseq"])
