#!/bin/bash
for i in $(eval echo {$1..$2})
  do 
     ./$i template2d_z4.msh 0 > $i.out &
 done
