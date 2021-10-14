#!/bin/bash

root -l -b <<EOF
.L vnlinearfit_Pt_TnPDw.C 
.q
EOF

echo entering to the loop...

for pt in '3,4.5' '3,6.5' '4.5,6.5'
do
  root -l -b -q 'vnlinearfit_Pt_TnPDw.C('$pt',1.6,2.4,20,120)'
done

for pt in '6.5,7.5' '6.5,9' '7.5,9' '9,12' '12,15' '15,50'
do
  root -l -b -q 'vnlinearfit_Pt_TnPDw.C('$pt',0,2.4,20,120)'
done

for pt in '3,4.5' '3,6.5' '4.5,6.5'
do
  root -l -b -q 'vnlinearfit_Pt_TnPDw.C('$pt',1.6,2.4,40,80)'
done

for pt in '6.5,7.5' '6.5,9' '7.5,9' '9,12' '12,15' '15,50'
do
  root -l -b -q 'vnlinearfit_Pt_TnPDw.C('$pt',0,2.4,40,80)'
done

for cent in '0,20' '20,40' '40,60' '60,80' '80,100'
do
  root -l -b -q 'vnlinearfit_Cent_TnPDw.C(6.5,50,0,2.4,'$cent')'
done

for pt in '3,4.5' '3,6.5' '4.5,6.5'
do
  root -l -b -q 'vnlinearfit_Pt_TnPUp.C('$pt',1.6,2.4,20,120)'
done

for pt in '6.5,7.5' '6.5,9' '7.5,9' '9,12' '12,15' '15,50'
do
  root -l -b -q 'vnlinearfit_Pt_TnPUp.C('$pt',0,2.4,20,120)'
done

for pt in '3,4.5' '3,6.5' '4.5,6.5'
do
  root -l -b -q 'vnlinearfit_Pt_TnPUp.C('$pt',1.6,2.4,40,80)'
done

for pt in '6.5,7.5' '6.5,9' '7.5,9' '9,12' '12,15' '15,50'
do
  root -l -b -q 'vnlinearfit_Pt_TnPUp.C('$pt',0,2.4,40,80)'
done

for cent in '0,20' '20,40' '40,60' '60,80' '80,100'
do
  root -l -b -q 'vnlinearfit_Cent_TnPUp.C(6.5,50,0,2.4,'$cent')'
done
