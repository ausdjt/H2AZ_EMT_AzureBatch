# todo: have to update fragmentlength based on macs2 predictd results
# have to change filterOut bed file

computeGCBias -b $1 \
              --effectiveGenomeSize 1800000000 \
              -g ~/Data/RefGenomes/CanFam3.1/UCSC/canFam3.2bit \
              --fragmentLength $2 \
              --biasPlot $1.png \
              --GCbiasFrequenciesFile $1.txt \
              --filterOut ~/Data/Annotations/hg19/UCSC/GCbias_filter_out.bed \
