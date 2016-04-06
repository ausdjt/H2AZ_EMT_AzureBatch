#!/bin/bash

bamCoverage --bam $1 \
            --outFileName $2 \
            --outFileFormat bigwig \
            --normalizeUsingRPKM \
            --centerReads \
            --missingDataAsZero yes \
            --binSize 1 \
            --fragmentLength 147 \
            --ignoreDuplicates \
            --smoothLength 10
