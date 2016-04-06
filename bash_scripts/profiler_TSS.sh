#!/bin/bash

function callProfiler()
{
  eval matrixFile="$1"
  eval outFileName="$2"
  eval plotTitle="$3"
  eval yAxisLabel="$4"
  eval startLabel="$5"


  plotProfile  --matrixFile ${matrixFile} \
            --outFileName ${outFileName} \
            --plotType fill \
            --startLabel "${startLabel}" \
            --plotTitle "${plotTitle}" \
            --yAxisLabel "${yAxisLabel}" \
            --yMin 0 \
            --yMax 120
}

if [ $# -ne 5 ]; then
  echo $0: usage: profiler_TSS.sh matrixFile outFileName yAxisLabel
  exit 1
fi

var1=$1
var2=$2
var3=$3
var4=$4
var5=$5

callProfiler "\${var1}" "\${var2}" "\${var3}" "\${var4}" "\${var5}"
