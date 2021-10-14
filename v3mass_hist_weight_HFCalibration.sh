#!/bin/bash

root -l -b <<EOF
.L v3mass_hist_weight_HFDw_211012.C 
.q
EOF

echo entering to the loop...

#pt bin1 cent10-60
for pt in '3,4.5' '3,6.5' '4.5,6.5'
do
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C('$pt',1.6,2.4,20,120,-1)'
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C('$pt',1.6,2.4,20,120,0)'
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C('$pt',1.6,2.4,20,120,1)'
done

for pt in '6.5,9' '9,12' '12,50'
do
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C('$pt',0,2.4,20,120,-1)'
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C('$pt',0,2.4,20,120,0)'
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C('$pt',0,2.4,20,120,1)'
done

#pt bin2 cent20-40
for pt in '3,4.5' '3,6.5' '4.5,6.5'
do
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C('$pt',1.6,2.4,40,80,-1)'
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C('$pt',1.6,2.4,40,80,0)'
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C('$pt',1.6,2.4,40,80,1)'
done

for pt in '6.5,9' '9,12' '12,50'
do
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C('$pt',0,2.4,40,80,-1)'
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C('$pt',0,2.4,40,80,0)'
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C('$pt',0,2.4,40,80,1)'
done

#cent Bin
for cent in '0,40' '40,80' '80,180'
do
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C(6.5,50,0,2.4,'$cent',-1)'
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C(6.5,50,0,2.4,'$cent',0)'
  root -l -b -q 'v3mass_hist_weight_HFDw_211012.C(6.5,50,0,2.4,'$cent',1)'
done

#pt bin1 cent10-60
for pt in '3,4.5' '3,6.5' '4.5,6.5'
do
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C('$pt',1.6,2.4,20,120,-1)'
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C('$pt',1.6,2.4,20,120,0)'
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C('$pt',1.6,2.4,20,120,1)'
done

for pt in '6.5,9' '9,12' '12,50'
do
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C('$pt',0,2.4,20,120,-1)'
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C('$pt',0,2.4,20,120,0)'
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C('$pt',0,2.4,20,120,1)'
done

root -l -b <<EOF
.L v3mass_hist_weight_HFUp_211012.C 
.q
EOF

echo entering to the loop...

#pt bin2 cent20-40
for pt in '3,4.5' '3,6.5' '4.5,6.5'
do
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C('$pt',1.6,2.4,40,80,-1)'
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C('$pt',1.6,2.4,40,80,0)'
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C('$pt',1.6,2.4,40,80,1)'
done

for pt in '6.5,9' '9,12' '12,50'
do
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C('$pt',0,2.4,40,80,-1)'
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C('$pt',0,2.4,40,80,0)'
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C('$pt',0,2.4,40,80,1)'
done

#cent Bin
for cent in '0,40' '40,80' '80,180'
do
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C(6.5,50,0,2.4,'$cent',-1)'
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C(6.5,50,0,2.4,'$cent',0)'
  root -l -b -q 'v3mass_hist_weight_HFUp_211012.C(6.5,50,0,2.4,'$cent',1)'
done
