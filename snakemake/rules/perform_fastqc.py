__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-08-08"

from snakemake.exceptions import MissingInputException
import os

def getAllFASTQ(wildcards):
    fn = []
    for i in config[wildcards.assayID]:
        for j in config[wildcards.assayID][i]:
                fn.append("RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/fastq/" + j)
    return(fn)

rule dummy:
    input:
        "RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/processed_data/reports/"

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
        "{assayID}/{runID}/{processed_dir}/{reports_dir}/"
    shell:
        "/usr/local/bin/fastqc {input} --noextract --outdir  {output}"
