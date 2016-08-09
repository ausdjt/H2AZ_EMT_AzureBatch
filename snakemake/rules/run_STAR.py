  __author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-04-22"

from snakemake.exceptions import MissingInputException

wrapper_dir = "/home/skurscheid/Development/snakemake-wrappers/bio"

rule star_align_full:
    version:
        0.4
    params:
        runThreadN = config["STAR"]["runThreadN"]
    input:
        "./{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R1_001.QT.CA.fastq.gz",
        "./{assayID}/{runID}/{processed_dir}/trimmed_data/{unit}_R2_001.QT.CA.fastq.gz",
        index = lambda wildcards: config["references"]["STAR"][wildcards.reference_version]
    output:
        "./{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam"
    shell:
        """
            STAR --runMode alignReads \
                 --runThreadN {params.runThreadN} \
                 --genomeDir {input.index} \
                 --readFilesIn {input[0]} {input[1]} \
                 --readFilesCommand zcat \
                 --outTmpDir /home/skurscheid/tmp/{wildcards.unit} \
                 --outSAMmode Full \
                 --outSAMattributes Standard \
                 --outSAMtype BAM SortedByCoordinate \
                 --outStd BAM_SortedByCoordinate \
                 --alignEndsType EndToEnd\
                 > {output}
        """

rule bam_index_STAR_output:
    version:
        0.2
    input:
        "./{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam"
    output:
        "./{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam.bai"
    wrapper:
        "file://" + wrapper_dir + "/samtools/index/wrapper.py"

rule run_htseq_count:
    version:
        0.3
    params:
        htseq_dir = config["HTSeq_dir"],
        gtf = config["references"]["GTF"]
    input:
        bam = "./{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam",
        index = "./{assayID}/{runID}/{processed_dir}/{reference_version}/STAR/full/{unit}.aligned.bam.bai"
    output:
        "./{assayID}/{runID}/{processed_dir}/{reference_version}/HTSeq/count/{unit}.txt"
    shell:
        """
            {params.htseq_dir}/htseq-count --format=bam \
                                           --order=pos \
                                           --stranded=reverse \
                                           --type=exon \
                                           --idattr=gene_id \
                                           --order=pos \
                                           {input.bam} \
                                           {params.gtf} \
                                           > {output}
        """

# rule run_dexseq_count:
#     version:
#         0.1
#     params:
#         dexseq_dir = config["DEXSeq_dir"],
#         dex_gtf = config["references"]["DEX_GTF"]
#     input:
#         bam = "{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam",
#         index = "{outdir}/{reference_version}/STAR/full/{unit}.aligned.bam.bai"
#     output:
#         "{outdir}/{reference_version}/DEXSeq/count/{unit}.txt"
#     shell:
#         """
#             python {params.dexseq_dir}/dexseq_count.py --format=bam \
#                                                        --paired=yes \
#                                                        --order=pos \
#                                                        --stranded=reverse \
#                                                        {params.dex_gtf} \
#                                                        {input.bam} \
#                                                        {output}
#         """
#
#
# rule collect_insert_size_metrics:
#     version:
#         0.1
#     params:
#         sampling = config["Picard"]["sampling"]
#     input:
#         rules.star_align_full.output
#     output:
#         txt = "{outdir}/{reference_version}/PICARD/insert_size_metrics/{unit}.insert_size_metrics.txt",
#         pdf = "{outdir}/{reference_version}/PICARD/insert_size_metrics/{unit}.insert_size_metrics.pdf"
#     shell:
#         """
#             java -Djava.io.tmpdir=/home/skurscheid/tmp \
#             -Xmx36G \
#             -jar /home/skurscheid/Bioinformatics/picard-tools-1.131/picard.jar CollectInsertSizeMetrics \
#             I={input} \
#             O={output.txt} \
#             H={output.pdf} \
#             M=0.2
#         """
#
# rule run_rMats:
#     version:
#         0.1
#     params:
#         gtf = config["references"]["GTF"],
#         bin = "/home/skurscheid/Bioinformatics/rMATS.3.2.2.beta/RNASeq-MATS.py"
#     input:
#         getGroups
#     output:
#         "{outdir}/{reference_version}/rMATS/{tissue}/{condition}"
#     shell:
#         """
#             python {params.bin} -b1 {input[0]} \
#                                 -b2 {input[1]} \
#                                 -gtf {params.gtf} \
#                                 -t paired \
#                                 -len 76 \
#                                 -analysis U \
#                                 -libType fr-firststrand \
#                                 -o {output}
#         """
