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
WORKFLOWDIR = config["WORKFLOWDIR"]
wrapper_dir = home + WORKFLOWDIR + "snakemake-wrappers/bio"
include_prefix= home + WORKFLOWDIR + "H2AZ_EMT/snakemake/rules/"
ASSAY = config["ASSAY"] #RNA-Seq
RUNID = config["RUNID"] #NB501086_0082_RDomaschenz_JCSMR_mRNAseq

# includes for the actual scripts
include:
    include_prefix + "perform_cutadapt.py"
include:
    include_prefix + "run_kallisto.py"
include:
    include_prefix + "run_STAR.py"

rule all:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/kallisto/{unit}",
               assayID = ASSAY,
               runID = RUNID,
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"],
               unit = config["samples"][ASSAY][RUNID]),
        expand("{assayID}/{runID}/{outdir}/{reference_version}/HTSeq/count/{unit}.txt",
               assayID = ASSAY,
               runID = RUNID,
               outdir = config["processed_dir"],
               reference_version = config["references"]["CanFam3.1"]["version"],
               unit = config["samples"][ASSAY][RUNID])
