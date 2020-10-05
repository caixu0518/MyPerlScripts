#!/bin/bash

input=$1
output=${input}".zhuanzhi"


awk '{for(i=1;i<=NF;i++){a[FNR,i]=\$i}}END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]" "}print ""}}' ${input}  >  ${output}
