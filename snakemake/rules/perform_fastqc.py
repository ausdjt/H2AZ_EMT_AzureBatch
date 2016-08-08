__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-08-08"

from snakemake.exceptions import MissingInputException
import os

def getAllFASTQ(wildcards):
    fn = []
    for i in config[wildcards.assayID]:
            fn.append("./fastq/" + i)
    return(fn)

rule dummy:
    input:
        expand("{rdir}", rdir = config["reports_dir"], assayID = "RNA-Seq")

rule fastqc:
    version:
        0.2
    message:
        "Performing FastQC..."
    params:
        assayID = "RNA-Seq",
        rdir = config["reports_dir"],
        data_dir = config["raw_dir"]
    input:
        getAllFASTQ
    output:
        "reports"
    shell:
        "/usr/local/bin/fastqc {input} --noextract --outdir  {output}"
