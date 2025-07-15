#!bin/bash 

awk -F, '{for(i=1;i<=NF;i++){if($i ~ /pt:i:[0-9]+/){print $1, $i}}}' rows_with_pt.csv
