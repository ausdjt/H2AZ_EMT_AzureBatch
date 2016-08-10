__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-10"

from snakemake.exceptions import MissingInputException

wrapper_dir = "/home/skurscheid/Development/snakemake-wrappers/bio"

localrules:
    run_kallisto_uncompressed, run_kallisto_quant

rule run_kallisto_uncompressed:
    input:
        expand("./{assayID}/{runID}/{outdir}/{reference_version}/uncompressed/kallisto/{unit}",
               runID = "NB501086_0067_RDomaschenz_JCSMR_RNASeq",
               assayID = "RNA-Seq",
               outdir = config["processed_dir"],
               reference_version = config["references"]["version"],
               unit = config["RNA-Seq"])

rule run_kallisto_quant:
        input:
            expand("./{assayID}/{runID}/{outdir}/{reference_version}/kallisto/{unit}",
                   runID = "NB501086_0067_RDomaschenz_JCSMR_RNASeq",
                   assayID = "RNA-Seq",
                   outdir = config["processed_dir"],
                   reference_version = config["references"]["version"],
                   unit = config["RNA-Seq"])

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
        "./{assayID}/{runID}/{processed_dir}/{reference_version}/kallisto/{unit}"
    shell:
        """
            kallisto quant --index={input.ki} \
                           --output-dir={output} \
                           --threads=4 \
                           --bootstrap-samples={params.bootstraps} \
                           {input[0]} {input[1]}
        """

rule kallisto_quant_from_uncompressed:
    message:
        "Running kallisto quant with uncompressed FASTQ as input..."
    params:
        bootstraps = 1,
        threads = 4
    input:
        read1 = lambda wildcards: "./" + wildcards.assayID + "/" + wildcards.runID + "/fastq/" + config[wildcards.assayID][wildcards.unit][0].replace(".gz", ""),
        read2 = lambda wildcards: "./" + wildcards.assayID + "/" + wildcards.runID + "/fastq/" + config[wildcards.assayID][wildcards.unit][1].replace(".gz", ""),
        ki = lambda wildcards: config["references"]["kallisto"][wildcards.reference_version]
    output:
        "./{assayID}/{runID}/{processed_dir}/{reference_version}/uncompressed/kallisto/{unit}"
    shell:
        """
            kallisto quant --index={input.ki} \
                           --output-dir={output} \
                           --threads=4 \
                           --bootstrap-samples={params.bootstraps} \
                          {input.read1} {input.read2 }
        """
