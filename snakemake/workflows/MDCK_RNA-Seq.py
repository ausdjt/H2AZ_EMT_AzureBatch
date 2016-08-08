_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2015-04-22"

from snakemake.exceptions import MissingInputException

rule:
    version: 0.3

localrules:
    all

include_prefix="/home/skurscheid/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/snakemake/rules/"

include:
     include_prefix + "perform_fastqc.py"
include:
    include_prefix + "perform_cutadapt.py"
include:
    include_prefix + "run_kallisto.py"
include:
    include_prefix + "run_STAR.py"

rule all:
    input:
        expand("./{trim_data}/{unit}_{suffix}.QT.CA.fastq.gz",
               unit = config["RNA-Seq"],
               trim_data = config["trim_dir"],
               suffix = ["R1_001", "R2_001"]),
        expand("{outdir}/{reference_version}/kallisto/{unit}",
               outdir = config["processed_dir"],
               reference_version = config["references"]["version"],
               unit = config["RNA-Seq"]),
        expand("{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam",
               outdir = config["processed_dir"],
               reference_version = config["references"]["version"],
               unit = config["RNA-Seq"]),
        expand("{outdir}/{reference_version}/HTSeq/count/{unit}.txt",
               outdir = config["processed_dir"],
               reference_version = config["references"]["version"],
               unit = config["RNA-Seq"]),
