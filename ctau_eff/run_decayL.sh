#!/bin/bash

root -l -b <<EOF
.L psuedo_proper_decay_length.C()
.q
EOF

echo entering to the loop...

for PRNP in '0' '1'
do
  for y in '1.6,2.4' #'0,2,4'
  do
    for pt in '3,6.5' '6.5,8' '8,10' '10,15' '15,20' '20,30'
    do
      echo  $outputDir 'psuedo_proper_decay_length.C('20,120','$pt','$y','0','2.6,3.5','$PRNP')'
      root -l -q $outputDir 'psuedo_proper_decay_length.C('20,120','$pt','$y','0','2.6,3.5','4.0','$PRNP')'
    done
  done
done
