# run up to 16 gzips in parallel
parallel -j 16 "gzip {}" ::: *.fastq
