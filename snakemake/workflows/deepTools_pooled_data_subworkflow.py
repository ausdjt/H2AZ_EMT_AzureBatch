__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-12-05"

from snakemake.exceptions import MissingInputException
import snakemake.utils
import os
import pdb

rule:
    version: 0.1

# run parameters as variables
ASSAY = config["ASSAY"]
RUNID = config["RUNID"]
OUTDIR = config["processed_dir"]
REFVERSION = config["references"]["CanFam3.1"]["version"][0]
QUAL = config["alignment_quality"]
home = os.environ['HOME']
WORKFLOWDIR = config["WORKFLOWDIR"]
wrapper_dir = home + WORKFLOWDIR + "snakemake-wrappers/bio"
include_prefix= home + WORKFLOWDIR + "H2AZ_EMT/snakemake/rules/"
configfile: home + WORKFLOWDIR + "H2AZ_EMT/snakemake/configs/config.json"

# includes for the actual scripts
include:
    include_prefix + "run_deepTools_pooled_data.py"


# targets
rule run_plotProfile_pooled_replicates:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/allSamples_{plotType}.{mode}.{region}.{suffix}",
                assayID = ASSAY,
                runID = RUNID,
                outdir = OUTDIR,
                reference_version = REFVERSION,
                application = "deepTools",
                tool = "plotProfile",
                command = ["reference-point", "scale-regions"],
                duplicates = ["duplicates_marked", "duplicates_removed"],
                referencePoint = "TSS",
                plotType = "se",
                region = ["allGenes", "TanEMTup", "TanEMTdown", "qPCRGenesUp", "qPCRGenesDown", "random100up", "random100down"],
                mode = ["MNase", "normal"],
                suffix = ["pdf", "bed", "data"])
