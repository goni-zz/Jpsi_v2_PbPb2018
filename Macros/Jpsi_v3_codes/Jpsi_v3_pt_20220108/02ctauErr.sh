#!/bin/bash

root -l -b <<EOF
.L CtauErr_v2Bin.C 
.q
EOF

for v2 in '0' '0.3' '0.6' '0.9' '1.2' '1.5' '1.8' '2.1' '-0.3' '-0.6' '-0.9' '-1.2' '-1.5' '-1.8' '-2.1' '-2.7'
do
  root -l -b -q 'CtauErr_v2Bin.C('$v2',3.0,6.5,1.6,2.4,20,120)'
  root -l -b -q 'CtauErr_v2Bin.C('$v2',6.5,9.0,0.0,2.4,20,120)'
  root -l -b -q 'CtauErr_v2Bin.C('$v2',9.0,12.0,0.0,2.4,20,120)'
  root -l -b -q 'CtauErr_v2Bin.C('$v2',12.0,15.0,0.0,2.4,20,120)'
  root -l -b -q 'CtauErr_v2Bin.C('$v2',15.0,50.0,0.0,2.4,20,120)'
done
