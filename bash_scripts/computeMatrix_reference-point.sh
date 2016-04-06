#!/bin/bash

computeMatrix reference-point \
              --regionsFileName $1 \
              --scoreFileName $2 \
              --referencePoint TSS \
              --beforeRegionStartLength 1000 \
              --afterRegionStartLength 500 \
              --binSize 10 \
              --sortRegions no \
              --skipZeros \
              --outFileName $3 \
              --outFileNameMatrix ${3}.txt \
              --outFileSortedRegions ${3}_sorted_regions.bed
