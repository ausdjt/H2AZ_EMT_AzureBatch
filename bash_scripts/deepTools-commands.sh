#!/bin/bash

export data_dir=~/Data/Tremethick/EMT/initial_analysis/GenomeWide/H2AZ/processed_data/duplicates_marked
export script_dir=~/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/bash_scripts
export refPoints_dir=~/Data/Annotations/CanFam3/Ensembl
export regionFiles_dir=~/Data/Tremethick/EMT/initial_analysis/GenomeWide/deepTools/regionFiles
export chip_dir=~/Data/Tremethick/EMT/initial_analysis/GenomeWide/H2AZ/processed_data/duplicates_marked
export input_dir=~/Data/Tremethick/EMT/initial_analysis/GenomeWide/Input/processed_data/duplicates_marked
export bigwig_dir=~/Data/Tremethick/EMT/initial_analysis/GenomeWide/bigwig

# TGFb-treated
bigwigCompare --bigwig1 $chip_dir/H2AZ_TGFb.bw \
              --bigwig2 $input_dir/Input_TGFb.bw \
              --numberOfProcessors 6 \
              --ratio log2 \
              --outFileName $bigwig_dir/H2AZ_TGFb_Input_log2ratio.bw

bigwigCompare --bigwig1 $chip_dir/H2AZ_TGFb.bw \
              --bigwig2 $input_dir/Input_TGFb.bw \
              --numberOfProcessors 6 \
              --ratio subtract \
              --outFileName $bigwig_dir/H2AZ_TGFb_Input_subtracted.bw


bamCoverage --bam $chip_dir/H2AZ_TGFb.bam \
            --outFileName $bigwig_dir/H2AZ_TGFb_10bp_RPKM.bw \
            --MNase\
            --binSize 10 \
            --normalizeUsingRPKM\
            --smoothLength 30\
            --ignoreDuplicates\
            --centerReads

bamCoverage --bam $input_dir/Input_TGFb.bam \
            --outFileName $bigwig_dir/Input_TGFb_10bp_RPKM.bw \
            --MNase\
            --binSize 10 \
            --normalizeUsingRPKM\
            --smoothLength 30\
            --ignoreDuplicates\
            --centerReads

bamCoverage --bam $chip_dir/H2AZ_TGFb.bam \
            --outFileName $bigwig_dir/H2AZ_TGFb_1bp_RPKM.bw \
            --binSize 1 \
            --normalizeUsingRPKM\
            --smoothLength 10\
            --ignoreDuplicates\
            --centerReads

bamCoverage --bam $input_dir/Input_TGFb.bam \
            --outFileName $bigwig_dir/Input_TGFb_1bp_RPKM.bw \
            --binSize 1 \
            --normalizeUsingRPKM\
            --smoothLength 10\
            --ignoreDuplicates\
            --centerReads

${script_dir}/computeMatrix_reference-point.sh  ${refPoints_dir}/EMT_down_genes_TSS.bed \
                                                ${data_dir}/H2AZ_TGFb.bw \
                                                ${data_dir}/H2AZ_TGFb_EMT_down_genes_TSS

${script_dir}/profiler_TSS.sh                   ${data_dir}/H2AZ_TGFb_EMT_down_genes_TSS \
                                                ${data_dir}/H2AZ_TGFb_EMT_down_genes_TSS.pdf \
                                                "H2A.Z TGFb-treated EMT down" \
                                                "mean coverage"

${script_dir}/computeMatrix_reference-point.sh  ${refPoints_dir}/EMT_up_genes_TSS.bed \
                                                ${data_dir}/H2AZ_TGFb.bw \
                                                ${data_dir}/H2AZ_TGFb_EMT_up_genes_TSS

${script_dir}/profiler_TSS.sh                   ${data_dir}/H2AZ_TGFb_EMT_up_genes_TSS \
                                                ${data_dir}/H2AZ_TGFb_EMT_up_genes_TSS.pdf \
                                                "H2A.Z TGFb-treated EMT up" \
                                                "mean coverage" \
                                                "TSS"

# 2016-04-07 changed BAM files to use de-duplicated data
export script_dir=~/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/bash_scripts
export refPoints_dir=~/Data/References/Annotations/CanFam3/Ensembl
export chip_dir=~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_removed
export input_dir=~/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_removed
export bigwig_dir=~/Data/Tremethick/EMT/GenomeWide/bigwig

bamCoverage --bam $chip_dir/H2AZ_TGFb.bam \
            --outFileName $bigwig_dir/H2AZ_TGFb_10bp_RPKM.bw \
            --MNase\
            --binSize 10 \
            --normalizeUsingRPKM\
            --smoothLength 30\
            --ignoreDuplicates\
            --centerReads

bamCoverage --bam $input_dir/Input_TGFb.bam \
            --outFileName $bigwig_dir/Input_TGFb_10bp_RPKM.bw \
            --MNase\
            --binSize 10 \
            --normalizeUsingRPKM\
            --smoothLength 30\
            --ignoreDuplicates\
            --centerReads

bamCoverage --bam $chip_dir/H2AZ_TGFb.bam \
            --outFileName $bigwig_dir/H2AZ_TGFb_1bp_RPKM.bw \
            --binSize 1 \
            --normalizeUsingRPKM\
            --smoothLength 10\
            --ignoreDuplicates\
            --centerReads

bamCoverage --bam $input_dir/Input_TGFb.bam \
            --outFileName $bigwig_dir/Input_TGFb_1bp_RPKM.bw \
            --binSize 1 \
            --normalizeUsingRPKM\
            --smoothLength 10\
            --ignoreDuplicates\
            --centerReads




# TGFb-treated top 5 upregulated EMT genes
${script_dir}/computeMatrix_reference-point.sh  ${refPoints_dir}/top5_EMT_up.bed \
                                                ${data_dir}/H2AZ_TGFb.bw \
                                                ${data_dir}/top5_EMT_up

${script_dir}/profiler_TSS.sh                   ${data_dir}/top5_EMT_up \
                                                ${data_dir}/top5_EMT_up.pdf \
                                                "H2A.Z TGFb-treated EMT up Top5" \
                                                "mean coverage"


# WT
bigwigCompare --bigwig1 $chip_dir/H2AZ_WT.bw \
              --bigwig2 $input_dir/Input_WT.bw \
              --numberOfProcessors 6 \
              --ratio log2 \
              --outFileName $bigwig_dir/H2AZ_WT_Input_log2ratio.bw

bigwigCompare --bigwig1 $chip_dir/H2AZ_WT.bw \
              --bigwig2 $input_dir/Input_WT.bw \
              --numberOfProcessors 6 \
              --ratio subtract \
              --outFileName $bigwig_dir/H2AZ_WT_Input_subtracted.bw


bamCoverage --bam $chip_dir/H2AZ_WT.bam \
            --outFileName $bigwig_dir/H2AZ_WT_10bp_RPKM.bw \
            --MNase\
            --binSize 10 \
            --normalizeUsingRPKM\
            --smoothLength 30\
            --ignoreDuplicates\
            --centerReads

bamCoverage --bam $chip_dir/H2AZ_WT.bam \
            --outFileName $bigwig_dir/H2AZ_WT_10bp_RPKM.bw \
            --MNase\
            --binSize 10 \
            --normalizeUsingRPKM\
            --smoothLength 30\
            --ignoreDuplicates\
            --centerReads

bamCoverage --bam $input_dir/Input_WT.bam \
            --outFileName $bigwig_dir/Input_WT_10bp_RPKM.bw \
            --binSize 1 \
            --normalizeUsingRPKM\
            --smoothLength 10\
            --ignoreDuplicates\
            --centerReads

bamCoverage --bam $input_dir/Input_WT.bam \
            --outFileName $bigwig_dir/Input_WT_1bp_RPKM.bw \
            --binSize 1 \
            --normalizeUsingRPKM\
            --smoothLength 10\
            --ignoreDuplicates\
            --centerReads

${script_dir}/computeMatrix_reference-point.sh  ${refPoints_dir}/EMT_down_genes_TSS.bed \
                                                ${data_dir}/H2AZ_WT.bw \
                                                ${data_dir}/H2AZ_WT_EMT_down_genes_TSS

${script_dir}/profiler_TSS.sh                   ${data_dir}/H2AZ_WT_EMT_down_genes_TSS \
                                                ${data_dir}/H2AZ_WT_EMT_down_genes_TSS.pdf \
                                                "H2A.Z WT EMT down" \
                                                "mean coverage"

${script_dir}/computeMatrix_reference-point.sh  ${refPoints_dir}/EMT_up_genes_TSS.bed \
                                                ${data_dir}/H2AZ_WT.bw \
                                                ${data_dir}/H2AZ_WT_EMT_up_genes_TSS

${script_dir}/profiler_TSS.sh                   ${data_dir}/H2AZ_WT_EMT_up_genes_TSS \
                                                ${data_dir}/H2AZ_WT_EMT_up_genes_TSS.pdf \
                                                "H2A.Z WT EMT up" \
                                                "mean coverage"

${script_dir}/computeMatrix_reference-point.sh  ${refPoints_dir}/all_genes_TSS.bed \
                                                ${data_dir}/H2AZ_WT.bw \
                                                ${data_dir}/H2AZ_WT_all_genes_TSS

${script_dir}/computeMatrix_reference-point.sh  ${refPoints_dir}/all_genes_TSS.bed \
                                                ${data_dir}/H2AZ_TGFb.bw \
                                                ${data_dir}/H2AZ_TGFb_all_genes_TSS

${script_dir}/profiler_TSS.sh                   ${data_dir}/H2AZ_WT_all_genes_TSS \
                                                ${data_dir}/H2AZ_WT_all_genes_TSS.pdf \
                                                "H2A.Z WT all genes" \
                                                "mean coverage"

${script_dir}/profiler_TSS.sh                   ${data_dir}/H2AZ_TGFb_all_genes_TSS \
                                                ${data_dir}/H2AZ_TGFb_all_genes_TSS.pdf \
                                                "H2A.Z TGFb all genes" \
                                                "mean coverage"

# groups of genes based on qPCR data hierarchical clustering analysis
computeMatrix reference-point -S $bigwig_dir/H2AZ_TGFb_10bp_RPKM.bw \
                                 $bigwig_dir/H2AZ_WT_10bp_RPKM.bw \
                              -R $regionFiles_dir/cluster1.red.bed \
                                 $regionFiles_dir/cluster2.green.bed \
                                 $regionFiles_dir/cluster3.blue.bed \
                                 $regionFiles_dir/cluster4.yellow.bed \
                                 --beforeRegionStartLength 1000 \
                                 --afterRegionStartLength 1000 \
                                 --binSize 10 \
                                 --sortRegions no \
                                 --skipZeros \
                                 --outFileName clustered.gz

plotProfile -m clustered.gz \
            -out clusteredH2AZ.pdf \
            --plotType lines \
            --colors red green blue yellow \
            --perGroup \
            --numPlotsPerRow 1\

# groups of genes based on qPCR data hierarchical clustering analysis (SD >0.5)
computeMatrix reference-point -S $bigwig_dir/H2AZ_TGFb_10bp_RPKM.bw \
                                 $bigwig_dir/H2AZ_WT_10bp_RPKM.bw \
                              -R $regionFiles_dir/sd05.cluster1.red.bed \
                                 $regionFiles_dir/sd05.cluster2.green.bed \
                                 $regionFiles_dir/sd05.cluster3.blue.bed \
                                 --beforeRegionStartLength 1000 \
                                 --afterRegionStartLength 1000 \
                                 --binSize 10 \
                                 --sortRegions no \
                                 --skipZeros \
                                 --outFileName clustered_sd05.gz

plotProfile -m clustered_sd05.gz \
            -out clusteredH2AZ_sd05.pdf \
            --plotType lines \
            --colors red green blue yellow \
            --perGroup \
            --numPlotsPerRow 1\

# plotting ChIP-Input data
computeMatrix reference-point -S $bigwig_dir/H2AZ_TGFb_Input_subtracted.bw \
                                 $bigwig_dir/H2AZ_WT_Input_subtracted.bw \
                              -R $regionFiles_dir/sd05.cluster1.red.bed \
                                 $regionFiles_dir/sd05.cluster2.green.bed \
                                 $regionFiles_dir/sd05.cluster3.blue.bed \
                                 --beforeRegionStartLength 1000 \
                                 --afterRegionStartLength 1000 \
                                 --binSize 10 \
                                 --sortRegions no \
                                 --skipZeros \
                                 --outFileName clustered_sd05_ChIP-Input.gz

plotProfile -m clustered_sd05_ChIP-Input.gz \
            -out clusteredH2AZ_sd05_ChIP-Input.pdf \
            --plotType lines \
            --colors red green blue yellow \
            --perGroup \
            --numPlotsPerRow 1
