#!/bin/bash

#root -l -b <<EOF
#.L v2mass_hist_weight_nominal_210928_Up.C 
#.q
#EOF

echo entering to the loop...

#for ctaucut in '-1' '0' '1'
#do
#  root -l -b -q doSimultaneousV2MassFit_weight_pt12_15_cent20_120_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt15_50_cent20_120_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt3_4p5_cent20_120_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt3_6p5_cent20_120_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt4p5_6p5_cent20_120_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_7p5_cent20_120_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_9_cent20_120_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt7p5_9_cent20_120_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt9_12_cent20_120_Up.C'('$ctaucut')'
#done

for ctaucut in '-1' '0' '1'
do
  root -l -b -q doSimultaneousV2MassFit_weight_pt15_50_cent40_80_Up.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt12_50_cent40_80_Up.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt15_50_cent40_80_Up.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt3_4p5_cent40_80_Up.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt3_6p5_cent40_80_Up.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt4p5_6p5_cent40_80_Up.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_7p5_cent40_80_Up.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_9_cent40_80_Up.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt7p5_9_cent40_80_Up.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt9_12_cent40_80_Up.C'('$ctaucut')'
done

#for ctaucut in '-1' '0' '1'
#do
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_50_cent0_20_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_50_cent20_40_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_50_cent40_60_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_50_cent60_80_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_50_cent80_100_Up.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_50_cent100_180_Up.C'('$ctaucut')'
#done
#
#for ctaucut in '-1' '0' '1'
#do
#  root -l -b -q doSimultaneousV2MassFit_weight_pt12_15_cent20_120_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt15_50_cent20_120_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt3_4p5_cent20_120_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt3_6p5_cent20_120_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt4p5_6p5_cent20_120_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_7p5_cent20_120_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_9_cent20_120_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt7p5_9_cent20_120_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt9_12_cent20_120_Dw.C'('$ctaucut')'
#done

for ctaucut in '-1' '0' '1'
do
  root -l -b -q doSimultaneousV2MassFit_weight_pt15_50_cent40_80_Dw.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt12_50_cent40_80_Dw.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt15_50_cent40_80_Dw.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt3_4p5_cent40_80_Dw.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt3_6p5_cent40_80_Dw.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt4p5_6p5_cent40_80_Dw.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_7p5_cent40_80_Dw.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_9_cent40_80_Dw.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt7p5_9_cent40_80_Dw.C'('$ctaucut')'
  root -l -b -q doSimultaneousV2MassFit_weight_pt9_12_cent40_80_Dw.C'('$ctaucut')'
done

#for ctaucut in '-1' '0' '1'
#do
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_50_cent0_20_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_50_cent20_40_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_50_cent40_60_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_50_cent60_80_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_50_cent80_100_Dw.C'('$ctaucut')'
#  root -l -b -q doSimultaneousV2MassFit_weight_pt6p5_50_cent100_180_Dw.C'('$ctaucut')'
#done
