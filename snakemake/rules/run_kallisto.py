__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-10"

from snakemake.exceptions import MissingInputException

wrapper_dir = "/home/skurscheid/Development/snakemake-wrappers/bio"

rule kallisto_quant:
    message:
        "Running kallisto quant..."
    params:
        bootstraps = config["kallisto"]["bootstraps"],
        threads = 4
    input:
        rules.cutadapt_pe.output,
        ki = lambda wildcards: config["references"]["kallisto"][wildcards.reference_version]
    output:
        protected("./{assayID}/{runID}/{processed_dir}/{reference_version}/kallisto/{unit}")
    wrapper:
        "file://" + wrapper_dir + "/kallisto/quant/wrapper.py"
