__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-08-08"

from snakemake.exceptions import MissingInputException
import os

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
        lambda wildcards: config["samples"][wildcards.assayID][wildcards.runID][wildcards.sample]
    output:
        "{assayID}/{runID}/{processed_dir}/{reports_dir}/{sample}"
    shell:
        "/usr/local/bin/fastqc {input} --noextract --outdir {output}"
