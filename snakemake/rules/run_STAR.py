  __author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-22"

# this set of rules is meant to be imported by the master workflow document

from snakemake.exceptions import MissingInputException

rule star_align_full:
    version:
        0.4
    params:
        runThreadN = config["program_parameters"]["STAR"]["runThreadN"],
        trim_dir = config["trim_dir"],
        tmp = temp("{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/tmp/{unit}")
    input:
        read1 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R1_001.QT.CA.fastq.gz",
        read2 = "{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R2_001.QT.CA.fastq.gz",
        index = lambda wildcards: home + config["references"]["CanFam3.1"]["STAR"][wildcards.reference_version]
    output:
        bam = "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {params.runThreadN} \
                 --genomeDir {input.index} \
                 --readFilesIn {input.read1} {input.read2}\
                 --readFilesCommand zcat \
                 --outTmpDir {params.tmp} \
                 --outSAMmode Full \
                 --outSAMattributes Standard \
                 --outSAMtype BAM SortedByCoordinate \
                 --outStd BAM_SortedByCoordinate \
                 --alignEndsType EndToEnd\
                 > {output.bam}
        """

rule bam_index_STAR_output:
    version:
        0.2
    input:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam"
    output:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam.bai"
    wrapper:
        "file://" + wrapper_dir + "/samtools/index/wrapper.py"

rule run_htseq_count:
    version:
        0.3
    params:
        htseq_dir = config["HTSeq_dir"]
    input:
        bam = "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam",
        index = "{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam.bai",
        gtf = lambda wildcards: home + config["references"]["CanFam3.1"]["GTF"][wildcards.reference_version]
    output:
        "{assayID}/{runID}/{processed_dir}/{reference_version}/HTSeq/count/{unit}.txt"
    shell:
        """
            {params.htseq_dir}/htseq-count --format=bam \
                                           --order=pos \
                                           --stranded=reverse \
                                           --type=exon \
                                           --idattr=gene_id \
                                           --order=pos \
                                           {input.bam} \
                                           {input.gtf} \
                                           > {output}
        """
