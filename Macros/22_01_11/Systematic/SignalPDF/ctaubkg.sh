#!/bin/bash

root -l -b <<EOF
.L CtauBkg_v2Bin.C 
.q
EOF

for v2 in '0' '0.3' '0.6' '0.9' '1.2' '1.5' '1.8' '2.1' '2.4' '-0.3' '-0.6' '-0.9' '-1.2' '-1.5' '-1.8' '-2.1' '-2.4' '-2.7'
  do
  for cent in '0,20' '20,40' '40,60' '60,80' '80,100' '100,180'
  do
    root -l -b -q 'CtauBkg_v2Bin.C('$v2',6.5,50,0,2.4,'$cent')'
  done
done

for v2 in '0' '0.3' '0.6' '0.9' '1.2' '1.5' '1.8' '2.1' '2.4' '-0.3' '-0.6' '-0.9' '-1.2' '-1.5' '-1.8' '-2.1' '-2.4' '-2.7'
  do
  for pt in '3,4.5' '4.5,6.5' 3,6.5''
  do
    root -l -b -q 'CtauBkg_v2Bin.C('$v2','$pt',1.6,2.4,20,120)'
  done
done

for v2 in '0' '0.3' '0.6' '0.9' '1.2' '1.5' '1.8' '2.1' '2.4' '-0.3' '-0.6' '-0.9' '-1.2' '-1.5' '-1.8' '-2.1' '-2.4' '-2.7'
  do
  for pt in '6.5,7.5' '7.5,9' '9,12' '12,15' '15,50' '6.5,9'
  do
    root -l -b -q 'CtauBkg_v2Bin.C('$v2','$pt',0,2.4,20,120)'
  done
done
