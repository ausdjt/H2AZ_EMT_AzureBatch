__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-08-08"

from snakemake.exceptions import MissingInputException
import os

def getAllFASTQ(wildcards, assayID):
    fn = []
    for i in config[assayID]:
        for j in config[assayID][i]:
            fn.append("./fastq/" + i + "/" + j)
    return(fn)

rule fastqc:
    version:
        0.2
    message:
        "Performing FastQC..."
    params:
        rdir = config["reports_dir"],
        data_dir = config["raw_dir"]
    input:
        getAllFASTQ(assayID = "RNA-Seq")
    output:
        "reports"
    shell:
        "/usr/local/bin/fastqc {input} --noextract --outdir  {output}"
