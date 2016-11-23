STAR --runThreadN 64 \
     --runMode genomeGenerate \
     --genomeDir ./STAR_Index/ \
     --genomeFastaFiles ${HOME}/Data/References/Genomes/Canis_familiaris/Ensembl/Canis_familiaris.CanFam3.1.dna.toplevel.fa ~/Data/References/Transcriptomes/ERCC/ERCC92.fa \
     --sjdbGTFfile ${HOME}/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Canis_familiaris.CanFam3.1.84.gtf \
     --sjdbOverhang 76


STAR --runMode genomeGenerate \
     --runThreadN 30 \
     --genomeDir ./STAR_Index/ \
     --genomeFastaFiles ${HOME}/Data/References/Genomes/Canis_familiaris/Ensembl/Canis_familiaris.CanFam3.1.dna.toplevel.fa \
     --sjdbGTFfile ${HOME}/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Canis_familiaris.CanFam3.1.84.gtf \
     --sjdbOverhang 76
