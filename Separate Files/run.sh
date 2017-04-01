#!/bin/bash
for i in $(eval echo {0..$1})
  do 
     ./set$i template2d_z4.msh 0 > $i.out &
 done