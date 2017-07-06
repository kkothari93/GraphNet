#!/bin/bash
for i in $(eval echo {$1..$2})
  do 
     ./$i rect5.msh 0 > $i.out &
 done
