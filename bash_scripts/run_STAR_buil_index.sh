STAR --runThreadN 64 \
     --runMode genomeGenerate \
     --genomeDir ./STAR_Index/ \
     --genomeFastaFiles ~/Data/References/Genomes/Canis_familiaris/Ensembl/Canis_familiaris.CanFam3.1.dna.toplevel.fa ~/Data/References/Transcriptomes/ERCC/ERCC92.fa \
     --sjdbGTFfile ~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Canis_familiaris.CanFam3.1.84.gtf \
     --sjdbOverhang 76
