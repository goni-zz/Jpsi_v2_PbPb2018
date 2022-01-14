#!/bin/bash

root -l -b <<EOF
.L CtauBkg_v2Bin.C 
.q
EOF

for v2 in '0' '0.3' '0.6' '0.9' '1.2' '1.5' '1.8' '2.1' '2.4' '-0.3' '-0.6' '-0.9' '-1.2' '-1.5' '-1.8' '-2.1' '-2.4' '-2.7'
do
  root -l -b -q 'CtauBkg_v2Bin.C('$v2',6.5,50,0,2.4,100,180)'
done
