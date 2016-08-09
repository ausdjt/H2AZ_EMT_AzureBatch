i__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-02-27"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for trimming reads with cutadapt
(http://cutadapt.readthedocs.org/en/latest/guide.html#illumina-truseq)

For usage, include this in your workflow.
"""

rule cutadapt_pe:
    """Trims given paired-end reads with given parameters"""
    params:
        trim_params = config["trim_params"],
        trim_data = config["trim_dir"],
        raw_data = config["raw_dir"],
        cutadapt_dir = config["cutadapt_dir"]
    input:

    output:
        "./RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/processed_data/trimmed_data/{unit}_R1_001.QT.CA.fastq.gz",
        "./RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/processed_data/trimmed_data/{unit}_R2_001.QT.CA.fastq.gz"
    shell:
        """
            {params.cutadapt_dir}/cutadapt {params.trim_params} \
                                            -o {output[0]} \
                                            -p {output[1]} \
                                            {input[0]} \
                                            {input[1]}
        """

rule dummy_cutadapt:
     """Trim all reads with all supplied trimming parameters"""
     input:
         expand("./{trim_data}/{unit}_R1_001.QT.CA.fastq.gz", unit = config["RNA-Seq"], trim_data = config["trim_dir"]),
         expand("./{trim_data}/{unit}_R2_001.QT.CA.fastq.gz", unit = config["RNA-Seq"], trim_data = config["trim_dir"])
