#!/bin/bash
# 2015-12-11
# command requires samtools-0.1.7
cd /Users/u1001407/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT_147bp

# create subdirectories to make call to danpos2 cleaner
mkdir TGFb_H2AZ; mkdir TGFb_Input; mkdir WT_H2AZ; mkdir WT_Input

# link to source BAM files
# Input - pooled input experiments
ln -fs ~/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_marked/Input_WT_rep1_S5.Q10.sorted.MkDup.bam ./WT_Input/Input_WT_rep1.bam
ln -fs ~/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_marked/Input_WT_rep1_S5.Q10.sorted.MkDup.bam.bai ./WT_Input/Input_WT_rep1.bam.bai
ln -fs ~/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_marked/Input_WT_rep2_S6.Q10.sorted.MkDup.bam ./WT_Input/Input_WT_rep2.bam
ln -fs ~/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_marked/Input_WT_rep2_S6.Q10.sorted.MkDup.bam.bai ./WT_Input/Input_WT_rep2.bam.bai

ln -fs ~/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_marked/Input_TGFb_rep1_S7.Q10.sorted.MkDup.bam ./TGFb_Input/Input_TGFb_rep1.bam
ln -fs ~/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_marked/Input_TGFb_rep1_S7.Q10.sorted.MkDup.bam.bai ./TGFb_Input/Input_TGFb_rep1.bam.bai
ln -fs ~/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_marked/Input_TGFb_rep2_S8.Q10.sorted.MkDup.bam ./TGFb_Input/Input_TGFb_rep2.bam
ln -fs ~/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_marked/Input_TGFb_rep2_S8.Q10.sorted.MkDup.bam.bai ./TGFb_Input/Input_TGFb_rep2.bam.bai

# H2A.Z
# TGFb-treated
ln -s ~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/H2AZ_TGFb_rep1_S3.Q10.sorted.MkDup.bam ./TGFb_H2AZ/H2AZ_TGFb_rep1.bam
ln -s ~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/H2AZ_TGFb_rep1_S3.Q10.sorted.MkDup.bam.bai ./TGFb_H2AZ/H2AZ_TGFb_rep1.bam.bai
ln -s ~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/H2AZ_TGFb_rep2_S4.Q10.sorted.MkDup.bam ./TGFb_H2AZ/H2AZ_TGFb_rep2.bam
ln -s ~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/H2AZ_TGFb_rep2_S4.Q10.sorted.MkDup.bam.bai ./TGFb_H2AZ/H2AZ_TGFb_rep2.bam.bai

# H2A.Z
# WT


# danpos2 requires Python 2.7.x
source activate py27

# actual command
python ~/Bioinformatics/danpos-2.2.2/danpos.py dpos TGFb_H2AZ/:WT_H2AZ/ \
                                                    --bg TGFb_H2AZ/:TGFb_Input/,WT_H2AZ/:WT_Input/ \
                                                    --position_reference ~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT/reference_positions.xls \
                                                    --paired 1

# modified run
