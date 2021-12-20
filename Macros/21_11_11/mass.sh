#!/bin/bash

root -l -b <<EOF
.L MassFit_FixPar_Data_v2Bins.C 
.q
EOF

for v2 in '0' '0.3' '0.6' '0.9' '1.2' '1.5' '1.8' '2.1' '2.4' '-0.3' '-0.6' '-0.9' '-1.2' '-1.5' '-1.8' '-2.1' '-2.4' '-2.7'
do
  root -l -b -q 'MassFit_FixPar_Data_v2Bins.C('$v2',6.5,50,0,2.4,100,180)'
done
