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
        "./{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R1_001.QT.CA.fastq.gz",
        "./{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R2_001.QT.CA.fastq.gz",
        ki = lambda wildcards: config["references"]["kallisto"][wildcards.reference_version]
    output:
        protected("./{assayID}/{runID}/{processed_dir}/{reference_version}/kallisto/{unit}")
    shell:
        """
            kallisto quant --index={input.ki} \
                           --output-dir={output} \
                           --threads=4 \
                           --bootstrap-samples={params.bootstraps} \
                           {input[0]} {input[1]}
        """
