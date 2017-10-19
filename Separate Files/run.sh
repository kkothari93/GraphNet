#!/bin/bash
for ((i=$1; i<=$2; ++i));
  do 
     ./${i} $3 0 > ${i}.out &
 done
