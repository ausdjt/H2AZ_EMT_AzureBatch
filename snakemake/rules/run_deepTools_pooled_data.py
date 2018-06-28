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
    a = config["program_parameters"][wildcards["application"]]["computeMatrix"][wildcards["command"]]
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

rule computeMatrix_pooled_replicates_single_matrix:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
        program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        file = lambda wildcards: expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/bamCoverage/{mode}/{duplicates}/{sampleGroup}_{mode}_RPKM.bw",
                                       assayID = ASSAYID,
                                       runID = RUNID,
                                       outdir = OUTDIR,
                                       reference_version = REFVERSION,
                                       application = "deepTools",
                                       duplicates = wildcards["duplicates"],
                                       sampleGroup = ["H2AZ-WT", "H2AZ-TGFb", "Input-WT", "Input-TGFb"],
                                       region = wildcards["region"],
                                       mode = wildcards["mode"]),
        region = lambda wildcards: home + config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards["region"]]
    output:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/allSamples_{region}_{mode}.matrix.gz"
    shell:
        """
            {params.deepTools_dir}/computeMatrix {wildcards.command} \
                                                 --regionsFileName {input.region} \
                                                 --scoreFileName {input.file} \
                                                 --missingDataAsZero \
                                                 --skipZeros \
                                                 --numberOfProcessors {threads} \
                                                 {params.program_parameters} \
                                                 --outFileName {output.matrix_gz}
        """


rule plotProfile_pooled_replicates:
    version:
        0.1
    params:
        deepTools_dir = home + config["deepTools_dir"],
    input:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/allSamples_{region}_{mode}.matrix.gz"
    output:
        figure = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/allSamples_{plotType}.{mode}.{region}.pdf",
        data = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/allSamples_{plotType}.{mode}.{region}.data",
        regions = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/allSamples_{plotType}.{mode}.{region}.bed"
    shell:
        """
            {params.deepTools_dir}/plotProfile --matrixFile {input.matrix_gz} \
                                               --outFileName {output.figure} \
                                               --outFileNameData {output.data} \
                                               --outFileSortedRegions {output.regions} \
                                               --plotType {wildcards.plotType}
        """

rule computeMatrix_pooled_replicates:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
        program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        file = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/bamCoverage/{mode}/{duplicates}/merged_replicates/{sample_group}_{mode}_RPKM.bw",
        region = lambda wildcards: home + config["program_parameters"]["deepTools"]["regionFiles"][wildcards.reference_version][wildcards.region]
    output:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/{sample_group}_{region}_{mode}.matrix.gz"
    shell:
        """
            {params.deepTools_dir}/computeMatrix {wildcards.command} \
                                                 --regionsFileName {input.region} \
                                                 --scoreFileName {input.file} \
                                                 --missingDataAsZero \
                                                 --skipZeros \
                                                 --numberOfProcessors {threads} \
                                                 {params.program_parameters} \
                                                 --outFileName {output.matrix_gz}
        """

rule bam_coverage_pooled_replicates:
    version:
        0.1
    params:
        deepTools_dir = home + config["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        program_parameters = cli_parameters_bamCoverage
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        bam = "{assayID}/{runID}/{outdir}/{reference_version}/samtools/merge/{duplicates}/{sample_group}.bam"
    output:
        bigwig = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/bamCoverage/{mode}/{duplicates}/merged_replicates/{sample_group}_{mode}_{norm}.bw"
    shell:
        """
            {params.deepTools_dir}/bamCoverage --bam {input.bam} \
                                               --outFileName {output.bigwig} \
                                               --outFileFormat bigwig \
                                               {params.program_parameters} \
                                               --numberOfProcessors {threads} \
                                               --normalizeUsingRPKM \
                                               --ignoreForNormalization {params.ignore}
        """

rule bigwig_compare_pooled_replicates:
    version:
        0.1
    params:
        deepTools_dir = home + config["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"]
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        control = "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/bamCoverage/{mode}/{duplicates}/merged_replicates/{control}_{mode}_{norm}.bw",
        treatment = "{assayID}/{runID}/{outdir}/{reference_version}/deepTools/bamCoverage/{mode}/{duplicates}/merged_replicates/{treatment}_{mode}_{norm}.bw"
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{scaleFactors}/{treatment}_vs_{control}_{mode}_{ratio}_{norm}.bw"
    shell:
        """
            {params.deepTools_dir}/bigwigCompare --bigwig1 {input.treatment} \
                                                 --bigwig2 {input.control} \
                                                 --outFileName {output} \
                                                 --ratio {wildcards.ratio} \
                                                 --numberOfProcessors {threads}
        """
