#!/bin/bash
for i in $(eval echo {$1..$2})
  do 
     ./$i template2d.msh 0 > $i.out &
 done
