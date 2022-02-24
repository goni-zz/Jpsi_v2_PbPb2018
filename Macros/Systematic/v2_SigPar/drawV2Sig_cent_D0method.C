#include <iostream>
#include "../../Style.h"
#include "../../commonUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
//#include "../CMS_lumi_v2mass.h"
//#include "../CMS_lumi_square.C"
#include "../../tdrstyle.C"

using namespace std;

void drawV2Sig_cent_D0method(){

	setTDRStyle();
  gStyle->SetEndErrorSize(4);
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;

  const int nCentBins=6;
  float centBin[nCentBins+1]={0,10,20,30,40,50,90};
 
  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TGraphErrors* gPR = new TGraphErrors();
  TGraphErrors* gNP = new TGraphErrors();
  TGraphErrors* gPRLow = new TGraphErrors();
  TGraphErrors* gNPLow = new TGraphErrors();

  //TFile *fPRpt1 = new TFile("v2_pt_3_4.5.root","read");
  //TFile *fPRpt2 = new TFile("v2_pt_4.5_6.5.root","read");
  //TFile *fPRpt3 = new TFile("v2_pt_6.5_7.5.root","read");
  //TFile *fPRpt4 = new TFile("v2_pt_7.5_9.root","read");
  //TFile *fPRpt5 = new TFile("v2_pt_9_12.root","read");
  //TFile *fPRpt6 = new TFile("v2_pt_12_50.root","read");

  TFile *fCent1 = new TFile("v2_D0method_pt6.5-50.0_y0.0-2.4_muPt0.0_centrality0-20.root","read");
  TFile *fCent2 = new TFile("v2_D0method_pt6.5-50.0_y0.0-2.4_muPt0.0_centrality20-40.root","read");
  TFile *fCent3 = new TFile("v2_D0method_pt6.5-50.0_y0.0-2.4_muPt0.0_centrality40-60.root","read");
  TFile *fCent4 = new TFile("v2_D0method_pt6.5-50.0_y0.0-2.4_muPt0.0_centrality60-80.root","read");
  TFile *fCent5 = new TFile("v2_D0method_pt6.5-50.0_y0.0-2.4_muPt0.0_centrality80-100.root","read");
  TFile *fCent6 = new TFile("v2_D0method_pt6.5-50.0_y0.0-2.4_muPt0.0_centrality100-180.root","read");


  TH1D * PR1 = (TH1D*) fCent1 -> Get("hOutPR"); 
  TH1D * PR2 = (TH1D*) fCent2 -> Get("hOutPR"); 
  TH1D * PR3 = (TH1D*) fCent3 -> Get("hOutPR"); 
  TH1D * PR4 = (TH1D*) fCent4 -> Get("hOutPR"); 
  TH1D * PR5 = (TH1D*) fCent5 -> Get("hOutPR"); 
  TH1D * PR6 = (TH1D*) fCent6 -> Get("hOutPR"); 

  TH1D * NP1 = (TH1D*) fCent1 -> Get("hOutNP"); 
  TH1D * NP2 = (TH1D*) fCent2 -> Get("hOutNP"); 
  TH1D * NP3 = (TH1D*) fCent3 -> Get("hOutNP"); 
  TH1D * NP4 = (TH1D*) fCent4 -> Get("hOutNP"); 
  TH1D * NP5 = (TH1D*) fCent5 -> Get("hOutNP"); 
  TH1D * NP6 = (TH1D*) fCent6 -> Get("hOutNP"); 

  double xPRcent1, yPRcent1, xPRcent2, yPRcent2, xPRcent3, yPRcent3, xPRcent4, yPRcent4, xPRcent5, yPRcent5, xPRcent6, yPRcent6, xPRcent7, yPRcent7;
  double xPRcent1Err, yPRcent1Err, xPRcent2Err, yPRcent2Err, xPRcent3Err, yPRcent3Err, xPRcent4Err, yPRcent4Err, xPRcent5Err, yPRcent5Err, xPRcent6Err, yPRcent6Err, xPRcent7Err, yPRcent7Err;
 
  yPRcent1= PR1->GetBinContent(0);
  yPRcent2= PR2->GetBinContent(0);
  yPRcent3= PR3->GetBinContent(0);
  yPRcent4= PR4->GetBinContent(0);
  yPRcent5= PR5->GetBinContent(0);
  yPRcent6= PR6->GetBinContent(0);
  yPRcent1Err= PR1->GetBinError(0);
  yPRcent2Err= PR2->GetBinError(0);
  yPRcent3Err= PR3->GetBinError(0);
  yPRcent4Err= PR4->GetBinError(0);
  yPRcent5Err= PR5->GetBinError(0);
  yPRcent6Err= PR6->GetBinError(0);

  gPR->SetPoint(0,(centBin[0]+centBin[1])/2,yPRcent1);
  gPR->SetPointError(0,0,yPRcent1Err);
  gPR->SetPoint(1,(centBin[1]+centBin[2])/2,yPRcent2);
  gPR->SetPointError(1,0,yPRcent2Err);
  gPR->SetPoint(2,(centBin[2]+centBin[3])/2,yPRcent3);
  gPR->SetPointError(2,0,yPRcent3Err);
  gPR->SetPoint(3,(centBin[3]+centBin[4])/2,yPRcent4);
  gPR->SetPointError(3,0,yPRcent4Err);
  gPR->SetPoint(4,(centBin[4]+centBin[5])/2,yPRcent5);
  gPR->SetPointError(4,0,yPRcent5Err);
  gPR->SetPoint(5,(centBin[5]+centBin[6])/2,yPRcent6);
  gPR->SetPointError(5,0,yPRcent6Err);

  double xNPcent1, yNPcent1, xNPcent2, yNPcent2, xNPcent3, yNPcent3, xNPcent4, yNPcent4, xNPcent5, yNPcent5, xNPcent6, yNPcent6, xNPcent7, yNPcent7, yNPcent8;
  double xNPcent1Err, yNPcent1Err, xNPcent2Err, yNPcent2Err, xNPcent3Err, yNPcent3Err, xNPcent4Err, yNPcent4Err, xNPcent5Err, yNPcent5Err, xNPcent6Err, yNPcent6Err, xNPcent7Err, yNPcent7Err, yNPcent8Err;
 
  yNPcent1= NP1->GetBinContent(0);
  yNPcent2= NP2->GetBinContent(0);
  yNPcent3= NP3->GetBinContent(0);
  yNPcent4= NP4->GetBinContent(0);
  yNPcent5= NP5->GetBinContent(0);
  yNPcent6= NP6->GetBinContent(0);
  yNPcent1Err= NP1->GetBinError(0);
  yNPcent2Err= NP2->GetBinError(0);
  yNPcent3Err= NP3->GetBinError(0);
  yNPcent4Err= NP4->GetBinError(0);
  yNPcent5Err= NP5->GetBinError(0);
  yNPcent6Err= NP6->GetBinError(0);
//  
  //gNPLow->SetPoint(0,(centBin[0]+centBin[1])/2,yNPcent1);
  //gNPLow->SetPointError(0,(centBin[1]-centBin[0])/2,yNPcent1Err);
  //gNPLow->SetPoint(1,(centBin[1]+centBin[2])/2,yNPcent2);
  //gNPLow->SetPointError(1,(centBin[2]-centBin[1])/2,yNPcent2Err);
  //gNP->SetPoint(0,(centBin[2]+centBin[3])/2,yNPcent3);
  //gNP->SetPointError(0,(centBin[3]-centBin[2])/2,yNPcent3Err);
  //gNP->SetPoint(1,(centBin[3]+centBin[4])/2,yNPcent4);
  //gNP->SetPointError(1,(centBin[4]-centBin[3])/2,yNPcent4Err);
  //gNP->SetPoint(2,(centBin[4]+centBin[5])/2,yNPcent5);
  //gNP->SetPointError(2,(centBin[5]-centBin[4])/2,yNPcent5Err);
  //gNP->SetPoint(3,(centBin[5]+centBin[6])/2,yNPcent6);
  //gNP->SetPointError(3,(centBin[6]-centBin[5])/2,yNPcent6Err);

  //gNPLow->SetPoint(0,(ptBin[0]+ptBin[1])/2,yNPpt1);
  //gNPLow->SetPointError(0,0,yNPpt1Err);
  //gNPLow->SetPoint(1,(ptBin[1]+ptBin[2])/2,yNPpt2);
  //gNPLow->SetPointError(1,0,yNPpt2Err);
  //cout << "ERROR : " << yNPpt2Err << endl;
  //gNP->SetPoint(0,(ptBin[2]+ptBin[3])/2,yNPpt3);
  //gNP->SetPointError(0,0,yNPpt3Err);
  //gNP->SetPoint(1,(ptBin[3]+ptBin[4])/2,yNPpt4);
  //gNP->SetPointError(1,0,yNPpt4Err);
  //gNP->SetPoint(2,(ptBin[4]+ptBin[5])/2,yNPpt5);
  //gNP->SetPointError(2,0,yNPpt5Err);
  //gNP->SetPoint(3,(ptBin[5]+ptBin[6])/2,yNPpt6);
  //gNP->SetPointError(3,0,yNPpt6Err);
  //gNP->SetPoint(4,(ptBin[6]+ptBin[7])/2,yNPpt7);
  //gNP->SetPointError(4,0,yNPpt7Err);

  gNP->SetPoint(0,(centBin[0]+centBin[1])/2,yNPcent1);
  gNP->SetPointError(0,0,yNPcent1Err);
  gNP->SetPoint(1,(centBin[1]+centBin[2])/2,yNPcent2);
  gNP->SetPointError(1,0,yNPcent2Err);
  gNP->SetPoint(2,(centBin[2]+centBin[3])/2,yNPcent3);
  gNP->SetPointError(2,0,yNPcent3Err);
  gNP->SetPoint(3,(centBin[3]+centBin[4])/2,yNPcent4);
  gNP->SetPointError(3,0,yNPcent4Err);
  gNP->SetPoint(4,(centBin[4]+centBin[5])/2,yNPcent5);
  gNP->SetPointError(4,0,yNPcent5Err);
  gNP->SetPoint(5,(centBin[5]+centBin[6])/2,yNPcent6);
  gNP->SetPointError(5,0,yNPcent6Err);

  gPR->GetXaxis()->SetTitle("Centrality (%)");
  gPR->GetXaxis()->CenterTitle();
  gPR->GetYaxis()->SetTitle("v_{2}");
  gPR->GetYaxis()->CenterTitle();
  gPR->GetXaxis()->SetLimits(0.,100);
  gPR->SetMinimum(-0.08);
  gPR->SetMaximum(0.25);

  gPR->SetLineColor(kMagenta+2);
  gPR->SetMarkerColor(kMagenta+2);
  gPR->SetMarkerSize(2);
  gPRLow->SetLineColor(kBlue+2);
  gPRLow->SetMarkerColor(kBlue+2);
  gPRLow->SetMarkerSize(2);

  gNP->SetLineColor(kOrange+5);
  gNP->SetMarkerColor(kOrange+5);
  gNP->SetMarkerSize(3);
  gNP->SetMarkerStyle(33);

  gNPLow->SetLineColor(kRed+2);
  gNPLow->SetMarkerColor(kRed+2);
  gNPLow->SetMarkerSize(3);
  gNPLow->SetMarkerStyle(33);


  float pos_x = 0.12;
  float pos_x_mass = 0.25;
  float pos_y = 0.93;
  float pos_y_diff = 0.071;
  int text_color = 1;
  float text_size = 25;
  TString perc = "%";

  TLegend *leg = new TLegend(0.58,0.6,0.8,0.8);
  SetLegendStyle(leg);
  leg->AddEntry(gPR,"Prompt J/#psi, |y^{#mu#mu}| < 2.4","pl");
  leg->AddEntry(gNP,"Nonprompt J/#psi, |y^{#mu#mu}| < 2.4","pl");


  TCanvas *c1 = new TCanvas("c1","",900,800);
  c1->cd();
  gPR->Draw("APE");
  gNP->Draw("pe");
  gPRLow->Draw("pe");
  gNPLow->Draw("pe");
  leg->Draw("sames");
  jumSun(0,0,100,0);

  drawText("Scalar Product", pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  drawText("6.5 - 50 GeV/c", pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  //drawText("Cent. 10-60 %", pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);

  //leg = new TLegen5(0.74, 0.70, 0.925, 0.77);
  //SetLegendStyle(leg);
  //leg->SetTextSize(0.044);
  //leg->AddEntry(hist_v2,"J/#psi","pe");
  //leg->Draw("same");
  //globtex->DrawLatex(0.23, sz_init, "1.6 < |y| < 2.4");
  //globtex->DrawLatex(0.23, sz_init-sz_step, Form("Cent. 10-60%s",perc.Data()));
  CMS_lumi_v2mass(c1,iPeriod,iPos);
  c1->Update();
  c1->SaveAs(Form("new_v2_cent.pdf"));
  //c1->SaveAs("v2_pt_Cent1060.pdf");
  //TFile* fSys2S = new TFile("../Systematic/merged_sys.root","read");
  //TH1D* hPt1S = (TH1D*) fSys->Get("hSys_pt_merged");
  //TH1D* hPt2S = (TH1D*) fSys2S->Get("hSys_pt_merged");
  //
  //double yPt1Ssys = hPt1S->GetBinContent(1);
  //double yPt2Ssys = hPt2S->GetBinContent(1);

  //cout << "Ptegr1ated 1S : " << Form("%.3f",yPt1S) << " +/- " << Form("%.3f",yPt1SErr) << " +/- " << Form("%.3f",yPt1Ssys) << endl;
  //cout << "Ptegr1ated 2S : " << Form("%.3f",yPt2S) << " +/- " << Form("%.3f",yPt2SErr) << " +/- " << Form("%.3f",yPt2Ssys) << endl;

  TFile *out = new TFile(Form("out_new_v2_vs_cent.root"),"recreate");
//  out->cd();
//  final_v2_int_g->SetName("gr1_point_v2_vs_pt_Cent1060");
//  gsys_int->SetName("gr1_sys_v2_vs_pt_Cent1060");
//  final_v2_int_g->Write();
//  gsys_int->Write(); 
//  for(int i=0;i<nCentBin;i++){
//    gsys[i]->Write();
//    final_v2_g[i]->Write();
//  }
//	return;
} 
