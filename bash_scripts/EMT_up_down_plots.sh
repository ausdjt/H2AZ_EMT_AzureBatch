export bigwig_dir="/home/sebastian/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/deepTools/bamCompare/normal/duplicates_marked/SES"
export regionFiles_dir="/home/sebastian/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84"

computeMatrix reference-point --scoreFileName $bigwig_dir/H2AZ-TGFb_vs_Input-TGFb_normal_log2_RPKM.bw \
                              --regionsFileName $regionFiles_dir/allGenes.bed \
                                                $regionFiles_dir/Tan_et_al_EMT_down_genes.bed \
                                                $regionFiles_dir/Tan_et_al_EMT_down_genes.bed \
                              --beforeRegionStartLength 1500 \
                              --afterRegionStartLength 1500 \
                              --binSize 10 \
                              --sortRegions no \
                              --skipZeros \
                              --outFileName H2AZ-TGFb_vs_Input-TGFb_normal_log2_RPKM.gz

plotProfile -m H2AZ-TGFb_vs_Input-TGFb_normal_log2_RPKM.gz \
            -out H2AZ-TGFb_vs_Input-TGFb_normal_log2_RPKM.pdf \
            --plotType se \
            --colors red green blue yellow \
            --perGroup \
            --numPlotsPerRow 1\
