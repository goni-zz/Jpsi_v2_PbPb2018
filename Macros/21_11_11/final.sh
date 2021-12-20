#!/bin/bash

root -l -b <<EOF
.L Final2DFit_v2Bins.C 
.q
EOF

for v2 in '0' '0.3' '0.6' '0.9' '1.2' '1.5' '1.8' '2.1' '-0.3' '-0.6' '-0.9' '-1.2' '-1.5' '-1.8' '-2.1' '-2.4'
do
  root -l -b -q 'Final2DFit_v2Bins.C('$v2',6.5,9,0,2.4,20,120)'
done
