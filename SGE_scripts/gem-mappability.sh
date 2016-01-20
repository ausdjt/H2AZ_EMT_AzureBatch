#!/bin/bash
#$ -cwd
#$ -j y
#$ -V

export data_dir=~/Data/RefGenomes/Canis_familiaris/Ensembl

gem-mappability -I ${data_dir}/CanFam3.1.gem \
                -o ${data_dir}/CanFam3.1_100bp_mappability \
                -l 100 \
                -T 4
                
