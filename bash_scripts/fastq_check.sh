awk -F ":" '$1 ~ /^@/ {mid=$10;total++;count[mid]++}END{for(mid in count){if(!max||count[mid]>max) {max=count[mid];maxMid=mid};if(count[mid]==1){unique++}}; print total,unique,maxMid,count[maxMid]}'

cat Renae7_S7_R1_001.fastq | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'
