#!/bin/bash

root -l -b <<EOF
.L v3mass_hist_weight_HFDown.C 
.q
EOF

echo entering to the loop...

#pt='6.5,7.5'
#ctauLo='0.008'
#ctauHi='0.116'
#ctauLo='0.003'
#ctauHi='0.101'
#for pt in '3,4.5' '4.5,6.5'
#do
#  root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',1.6,2.4,40,80,-1)'
#  root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',1.6,2.4,40,80,0)'
#  root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',1.6,2.4,40,80,1)'
#done

# for pt in '3,4.5' '3,6.5' '4.5,6.5'
# do
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',1.6,2.4,20,120,-1)'
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',1.6,2.4,20,120,0)'
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',1.6,2.4,20,120,1)'
# done
#
# for pt in '6.5,7.5' '6.5,9' '7.5,9' '9,12' '12,15' '15,50'
# do
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',0,2.4,20,120,-1)'
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',0,2.4,20,120,0)'
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',0,2.4,20,120,1)'
# done
#
# for pt in '3,4.5' '3,6.5' '4.5,6.5'
# do
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',1.6,2.4,40,80,-1)'
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',1.6,2.4,40,80,0)'
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',1.6,2.4,40,80,1)'
# done
#
# for pt in '6.5,7.5' '6.5,9' '7.5,9' '9,12' '12,15' '15,50'
# do
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',0,2.4,40,80,-1)'
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',0,2.4,40,80,0)'
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C('$pt',0,2.4,40,80,1)'
# done
#
# for cent in '0,20' '20,40' '0,40' '40,60' '60,80' '80,100' '100,180'
# do
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C(6.5,50,0,2.4,'$cent',-1)'
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C(6.5,50,0,2.4,'$cent',0)'
  # root -l -b -q 'v3mass_hist_weight_AccSys_211013.C(6.5,50,0,2.4,'$cent',1)'
# done

for cent in '40,80'
do
  root -l -b -q 'v3mass_hist_weight_AccSys_211013.C(6.5,50,0,2.4,'$cent',-1)'
  root -l -b -q 'v3mass_hist_weight_AccSys_211013.C(6.5,50,0,2.4,'$cent',0)'
  root -l -b -q 'v3mass_hist_weight_AccSys_211013.C(6.5,50,0,2.4,'$cent',1)'
done
