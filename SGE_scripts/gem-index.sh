#!/bin/bash
#$ -cwd
#$ -j y
#$ -V

export data_dir=~/Data/RefGenomes/Canis_familiaris/Ensembl

gem-indexer -i ${data_dir}/Canis_familiaris.CanFam3.1.dna_rm.toplevel.fa \
            -o ${data_dir}/CanFam3.1 \
            
