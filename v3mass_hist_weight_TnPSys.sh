#!/bin/bash

root -l -b <<EOF
.L v3mass_hist_weight_TnPSysDw_210928.C
.q
EOF

echo entering to the loop...


for pt in '3,6.5'
do
 root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C('$pt',1.6,2.4,20,120,-1)'
 root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C('$pt',1.6,2.4,20,120,0)'
 root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C('$pt',1.6,2.4,20,120,1)'
done

for pt in '6.5,9' '9,12' '12,15' '15,50'
do
 root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C('$pt',0,2.4,20,120,-1)'
 root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C('$pt',0,2.4,20,120,0)'
 root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C('$pt',0,2.4,20,120,1)'
done

#for pt in '3,6.5'
#do
# root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C('$pt',1.6,2.4,40,80,-1)'
# root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C('$pt',1.6,2.4,40,80,0)'
# root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C('$pt',1.6,2.4,40,80,1)'
#done
#
#for pt in '6.5,9' '9,12' '12,15' '15,50'
##for pt in '12,15' '15,50'
#do
#  root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C('$pt',0,2.4,40,80,-1)'
#  root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C('$pt',0,2.4,40,80,0)'
#  root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C('$pt',0,2.4,40,80,1)'
#done

for cent in '0,40' '40,80'
do
  root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C(6.5,50,0,2.4,'$cent',-1)'
  root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C(6.5,50,0,2.4,'$cent',0)'
  root -l -b -q 'v3mass_hist_weight_TnPSysDw_210928.C(6.5,50,0,2.4,'$cent',1)'
done

root -l -b <<EOF
.L v3mass_hist_weight_TnPSysUp_210928.C
.q
EOF

for pt in '3,6.5'
do
 root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C('$pt',1.6,2.4,20,120,-1)'
 root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C('$pt',1.6,2.4,20,120,0)'
 root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C('$pt',1.6,2.4,20,120,1)'
done

for pt in '6.5,9' '9,12' '12,15' '15,50'
do
 root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C('$pt',0,2.4,20,120,-1)'
 root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C('$pt',0,2.4,20,120,0)'
 root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C('$pt',0,2.4,20,120,1)'
done

#for pt in '3,6.5'
#do
# root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C('$pt',1.6,2.4,40,80,-1)'
# root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C('$pt',1.6,2.4,40,80,0)'
# root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C('$pt',1.6,2.4,40,80,1)'
#done
#
#for pt in '6.5,9' '9,12' '12,15' '15,50'
##for pt in '12,15' '15,50'
#do
#  root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C('$pt',0,2.4,40,80,-1)'
#  root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C('$pt',0,2.4,40,80,0)'
#  root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C('$pt',0,2.4,40,80,1)'
#done

for cent in '0,40' '40,80'
do
  root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C(6.5,50,0,2.4,'$cent',-1)'
  root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C(6.5,50,0,2.4,'$cent',0)'
  root -l -b -q 'v3mass_hist_weight_TnPSysUp_210928.C(6.5,50,0,2.4,'$cent',1)'
done
