#!/bin/bash

root -l -b <<EOF
.L v2mass_hist_weight.C 
.q
EOF

echo entering to the loop...

for pt in '6.5,7.5'
do
  for ctau in '0' '0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1'
  do
    root -l -b -q 'v2mass_hist_weight.C('$pt',0,2.4,20,120,'$ctau',-1)'
    root -l -b -q 'v2mass_hist_weight.C('$pt',0,2.4,20,120,'$ctau',1)'
  done
done
