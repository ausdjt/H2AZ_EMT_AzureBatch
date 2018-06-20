__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-01-17"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os

"""
Rules for running deepTools analysis on ChIP-Seq data
For usage, include this in your workflow.
"""

def cli_parameters_computeMatrix(wildcards):
    a = config["program_parameters"][wildcards["application"]]["computeMatrix"]][wildcards["command"]]
    if wildcards["command"] == "reference-point":
        a["--referencePoint"] = wildcards.referencePoint
    return(a)

def cli_parameters_bamCoverage(wildcards):
    a = config["program_parameters"][wildcards["application"]]["bamCoverage"][wildcards["mode"]]
    b = str()
    for (key, val) in a.items():
        if val == " ":
            f = key + " "
            b = b + f
        else:
            f = key + "=" + val + " "
            b = b + f
    if wildcards["mode"] == "MNase":
        b = b + "--MNase"
    return(b.rstrip())

def get_computeMatrix_input(wildcards):
    fn = []
    path = "/".join((wildcards["assayID"],
                     wildcards["runID"],
                     config["processed_dir"],
                     config["references"]["CanFam3.1"]["version"][0],
                     wildcards["application"],
                     "bamCoverage",
                     wildcards["mode"],
                     wildcards["duplicates"]))
    for i in config["samples"][wildcards["assayID"]][wildcards["runID"]]:
        fn.append("/".join((path, "_".join((i, wildcards["mode"], "RPKM.bw")))))
    return(fn)

rule bamCoverage:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        program_parameters = cli_parameters_bamCoverage
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        bam = lambda wildcards: expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{duplicates}/{sample}.Q{qual}.{suffix}",
                                       assayID = wildcards["assayID"],
                                       runID = wildcards["runID"],
                                       outdir = wildcards["outdir"],
                                       reference_version = wildcards["reference_version"],
                                       application = "bowtie2",
                                       duplicates = wildcards["duplicates"],
                                       sample = wildcards["sample"],
                                       qual = config["alignment_quality"],
                                       suffix = "sorted.bam")
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{sample}_{mode}_{norm}.bw"
    shell:
        """
            {params.deepTools_dir}/bamCoverage --bam {input.bam} \
                                               --outFileName {output} \
                                               --outFileFormat bigwig \
                                               {params.program_parameters} \
                                               --numberOfProcessors {threads} \
                                               --normalizeUsingRPKM \
                                               --ignoreForNormalization {params.ignore}
        """

rule computeMatrix:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
        program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        file = get_computeMatrix_input,
        region = lambda wildcards: home + config["program_parameters"]["deepTools"]["regionFiles"][wildcards.region]
    output:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{region}_{mode}.matrix.gz"
    wrapper:
        "file://" + wrapper_dir + "/deepTools/computeMatrix/wrapper.py"


rule plotProfile:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
    input:
        rules.computeMatrix.output
    output:
        figure = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{region}.pdf",
        data = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{region}.data",
        regions = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{region}.bed"
    shell:
        """
            {params.deepTools_dir}/plotProfile --matrixFile {input.matrix_gz} \
                                               --outFileName {output.figure} \
                                               --outFileNameData {output.data} \
                                               --outFileSortedRegions {output.regions} \
                                               --plotType {wildcards.plotType}
        """
