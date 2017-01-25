import json
from pprint import pprint

with open("/Users/u1001407/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/snakemake/configs/config.json") as data_file:
    config = json.load(data_file)

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


wildcards = {"assayID" : "ChIP-Seq", "runID" : "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq", "application" : "deepTools", "mode" : "normal", "duplicates" : "duplicates_marked"}

wildcards = {"outdir" : "processed_data",
             "application" : "deepTools",
             "reference_version" : "CanFam3.1_ensembl84_ERCC",
             "runID" : "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
             "tool" : "computeMatrix",
             "duplicates" : "duplicates_removed",
             "command" : "scale-regions",
             "region" : "allGenes",
             "assayID" : "ChIP-Seq",
             "mode" : "normal",
             "referencePoint" : "TSS"}

wildcards = {"outdir" : "processed_data",
             "application" : "deepTools",
             "reference_version" : "CanFam3.1_ensembl84_ERCC",
             "runID" : "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
             "tool" : "computeMatrix",
             "duplicates" : "duplicates_removed",
             "command" : "scale-regions",
             "region" : "allGenes",
             "assayID" : "ChIP-Seq",
             "mode" : "MNase",
             "referencePoint" : "TSS"}

get_computeMatrix_input(wildcards)


file = expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{sample}_{mode}_{norm}.bw",
                      assayID = "ChIP-Seq",
                      runID = "NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq",
                      outdir = config["processed_dir"],
                      reference_version = config["references"]["CanFam3.1"]["version"][0],
                      application = "deepTools",
                      tool = "bamCoverage",
                      mode = ["MNase", "normal"],
                      duplicates = ["duplicates_marked", "duplicates_removed"],
                      sample = config["samples"]["ChIP-Seq"]["NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq"],
                      norm = "RPKM"),
b = "/".join("test", "two")
