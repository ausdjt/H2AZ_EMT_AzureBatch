__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-10"

# this set of rules is meant to be imported by the master workflow document

from snakemake.exceptions import MissingInputException

rule kallisto_quant:
    params:
        bootstraps = config["program_parameters"]["kallisto"]["bootstraps"],
        threads = 4,
        trim_dir = config["trim_dir"]
    input:
        read1 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R1_001.QT.CA.fastq.gz",
        read2 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R2_001.QT.CA.fastq.gz",
        ki = lambda wildcards: home + config["references"]["CanFam3.1"]["kallisto"][wildcards.reference_version]
    output:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/kallisto/{unit}"
    shell:
        """
            kallisto quant --index={input.ki} \
                           --output-dir={output} \
                           --threads=4 \
                           --bootstrap-samples={params.bootstraps} \
                           {input.read1} {input.read2}
        """
