from snakemake.exceptions import MissingInputException
import os

home = os.environ['HOME']

# run parameters as variables
RUNID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"
ASSAYID = "ChIP-Seq"
OUTDIR = config["processed_dir"]
REFVERSION = config["references"]["CanFam3.1"]["version"][0]
QUALITY = config["alignment_quality"]

def cli_parameters_computeMatrix(wildcards):
    a = config["program_parameters"][wildcards["application"]]["computeMatrix"][wildcards["command"]]
    if wildcards["command"] == "reference-point":
        a["--referencePoint"] = wildcards.referencePoint
    return(a)

# pseudo rules for build targets
rule run_computeMatrix_pooled_replicates:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/{sampleGroup}_{region}_{mode}.matrix.gz",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               application = "deepTools",
               tool = "computeMatrix",
               command = ["reference-point", "scale-regions"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               referencePoint = "TSS",
               sampleGroup = ["H2AZ-WT", "H2AZ-TGFb", "Input-WT", "Input-TGFb"],
               region = ["allGenes", "TanEMTup", "TanEMTdown", "qPCRGenesUp", "qPCRGenesDown", "random100up", "random100down"],
               mode = ["MNase", "normal"])


rule run_computeMatrix_pooled_replicates_bigwigCompare_single_matrix_WT:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/bigwigCompare/{duplicates}/{referencePoint}/{treatment}_vs_{control}_normal.{scaleFactors}.{ratio}_{norm}_{region}_{mode}.matrix.gz",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               application = "deepTools",
               command = ["reference-point", "scale-regions"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               referencePoint = "TSS",
               treatment = "H2AZ-WT",
               control = "Input-WT",
               scaleFactors = ["SES", "readCount"],
               ratio = ["log2", "subtract"],
               norm = "RPKM",
               region = ["allGenes", "TanEMTup", "TanEMTdown", "qPCRGenesUp", "qPCRGenesDown", "random100up", "random100down"],
               mode = ["MNase", "normal"])

rule run_computeMatrix_pooled_replicates_bigwigCompare_single_matrix_TGFb:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/bigwigCompare/{duplicates}/{referencePoint}/{treatment}_vs_{control}_normal.{scaleFactors}.{ratio}_{norm}_{region}_{mode}.matrix.gz",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               application = "deepTools",
               command = ["reference-point"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               referencePoint = "TSS",
               treatment = "H2AZ-TGFb",
               control = "Input-TGFb",
               scaleFactors = ["readCount"],
               ratio = ["subtract"],
               norm = "RPKM",
               region = ["allGenes", "TanEMTup", "TanEMTdown", "qPCRGenesUp", "qPCRGenesDown", "random100up", "random100down"],
               mode = ["MNase", "normal"])


# Actual run rules
# rule computeMatrix_pooled_replicates:
#     version:
#         0.2
#     params:
#         deepTools_dir = home + config["deepTools_dir"],
#         program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
#     threads:
#         lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
#     input:
#         file = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/bamCoverage/{mode}/{duplicates}/{sampleGroup}_{mode}_RPKM.bw",
#         region = lambda wildcards: home + config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards["region"]]
#     output:
#         matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/{sampleGroup}_{region}_{mode}.matrix.gz"
#     shell:
#         """
#             {params.deepTools_dir}/computeMatrix {wildcards.command} \
#                                                  --regionsFileName {input.region} \
#                                                  --scoreFileName {input.file} \
#                                                  --missingDataAsZero \
#                                                  --skipZeros \
#                                                  --numberOfProcessors {threads} \
#                                                  {params.program_parameters} \
#                                                  --outFileName {output.matrix_gz}
#         """

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

rule computeMatrix_pooled_replicates_bigwigCompare_single_matrix:
    version:
        0.1
    params:
        deepTools_dir = home + config["deepTools_dir"],
        program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        file = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/bigwigCompare/normal/{duplicates}/{scaleFactors}/{treatment}_vs_{control}_normal_{ratio}_{norm}.bw",
        region = lambda wildcards: home + config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards["region"]]
    output:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/bigwigCompare/{duplicates}/{referencePoint}/{treatment}_vs_{control}_normal.{scaleFactors}.{ratio}_{norm}_{region}_{mode}.matrix.gz"
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
