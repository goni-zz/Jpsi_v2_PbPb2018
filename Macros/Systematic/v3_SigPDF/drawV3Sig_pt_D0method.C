#include <iostream>
#include "../../Style.h"
#include "../../commonUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
//#include "../CMS_lumi_v2mass.h"
//#include "../CMS_lumi_square.C"
#include "../../tdrstyle.C"

using namespace std;

void drawV3Sig_pt_D0method(){

  setTDRStyle();
  gStyle->SetEndErrorSize(4);
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;

  const int nPtBins=5;
  float ptBin[nPtBins+1]={3,6.5,9,12,15,50};
 
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

  TFile *fPt1 = new TFile("v2_D0method_pt3.0-6.5_y1.6-2.4_muPt0.0_centrality20-120.root","read");
  TFile *fPt2 = new TFile("v2_D0method_pt6.5-9.0_y0.0-2.4_muPt0.0_centrality20-120.root","read");
  TFile *fPt3 = new TFile("v2_D0method_pt9.0-12.0_y0.0-2.4_muPt0.0_centrality20-120.root","read");
  TFile *fPt4 = new TFile("v2_D0method_pt12.0-15.0_y0.0-2.4_muPt0.0_centrality20-120.root","read");
  TFile *fPt5 = new TFile("v2_D0method_pt15.0-50.0_y0.0-2.4_muPt0.0_centrality20-120.root","read");


  TH1D * PR1 = (TH1D*) fPt1 -> Get("hOutPR"); 
  TH1D * PR2 = (TH1D*) fPt2 -> Get("hOutPR"); 
  TH1D * PR3 = (TH1D*) fPt3 -> Get("hOutPR"); 
  TH1D * PR4 = (TH1D*) fPt4 -> Get("hOutPR"); 
  TH1D * PR5 = (TH1D*) fPt5 -> Get("hOutPR"); 

  TH1D * NP1 = (TH1D*) fPt1 -> Get("hOutNP"); 
  TH1D * NP2 = (TH1D*) fPt2 -> Get("hOutNP"); 
  TH1D * NP3 = (TH1D*) fPt3 -> Get("hOutNP"); 
  TH1D * NP4 = (TH1D*) fPt4 -> Get("hOutNP"); 
  TH1D * NP5 = (TH1D*) fPt5 -> Get("hOutNP"); 

  double xPRpt1, yPRpt1, xPRpt2, yPRpt2, xPRpt3, yPRpt3, xPRpt4, yPRpt4, xPRpt5, yPRpt5, xPRpt6, yPRpt6, xPRpt7, yPRpt7;
  double xPRpt1Err, yPRpt1Err, xPRpt2Err, yPRpt2Err, xPRpt3Err, yPRpt3Err, xPRpt4Err, yPRpt4Err, xPRpt5Err, yPRpt5Err, xPRpt6Err, yPRpt6Err, xPRpt7Err, yPRpt7Err;
 
  yPRpt1= PR1->GetBinContent(0);
  yPRpt2= PR2->GetBinContent(0);
  yPRpt3= PR3->GetBinContent(0);
  yPRpt4= PR4->GetBinContent(0);
  yPRpt5= PR5->GetBinContent(0);
  yPRpt1Err= PR1->GetBinError(0);
  yPRpt2Err= PR2->GetBinError(0);
  yPRpt3Err= PR3->GetBinError(0);
  yPRpt4Err= PR4->GetBinError(0);
  yPRpt5Err= PR5->GetBinError(0);
  cout << "15-50 Error : " << yPRpt5Err << endl;
//  
  //gPRLow->SetPoint(0,(ptBin[0]+ptBin[1])/2,yPRpt1);
  //gPRLow->SetPointError(0,(ptBin[1]-ptBin[0])/2,yPRpt1Err);
  //gPRLow->SetPoint(1,(ptBin[1]+ptBin[2])/2,yPRpt2);
  //gPRLow->SetPointError(1,(ptBin[2]-ptBin[1])/2,yPRpt2Err);
  //gPR->SetPoint(0,(ptBin[2]+ptBin[3])/2,yPRpt3);
  //gPR->SetPointError(0,(ptBin[3]-ptBin[2])/2,yPRpt3Err);
  //gPR->SetPoint(1,(ptBin[3]+ptBin[4])/2,yPRpt4);
  //gPR->SetPointError(1,(ptBin[4]-ptBin[3])/2,yPRpt4Err);
  //gPR->SetPoint(2,(ptBin[4]+ptBin[5])/2,yPRpt5);
  //gPR->SetPointError(2,(ptBin[5]-ptBin[4])/2,yPRpt5Err);
  //gPR->SetPoint(3,(ptBin[5]+ptBin[6])/2,yPRpt6);
  //gPR->SetPointError(3,(ptBin[6]-ptBin[5])/2,yPRpt6Err);
  gPRLow->SetPoint(0,(ptBin[0]+ptBin[1])/2,yPRpt1);
  gPRLow->SetPointError(0,0,yPRpt1Err);
  gPR->SetPoint(0,(ptBin[1]+ptBin[2])/2,yPRpt2);
  gPR->SetPointError(0,0,yPRpt2Err);
  gPR->SetPoint(1,(ptBin[2]+ptBin[3])/2,yPRpt3);
  gPR->SetPointError(1,0,yPRpt3Err);
  gPR->SetPoint(2,(ptBin[3]+ptBin[4])/2,yPRpt4);
  gPR->SetPointError(2,0,yPRpt4Err);
  gPR->SetPoint(3,(ptBin[4]+ptBin[5])/2,yPRpt5);
  gPR->SetPointError(3,0,yPRpt5Err);

  double xNPpt1, yNPpt1, xNPpt2, yNPpt2, xNPpt3, yNPpt3, xNPpt4, yNPpt4, xNPpt5, yNPpt5, xNPpt6, yNPpt6, xNPpt7, yNPpt7, yNPpt8;
  double xNPpt1Err, yNPpt1Err, xNPpt2Err, yNPpt2Err, xNPpt3Err, yNPpt3Err, xNPpt4Err, yNPpt4Err, xNPpt5Err, yNPpt5Err, xNPpt6Err, yNPpt6Err, xNPpt7Err, yNPpt7Err, yNPpt8Err;
 
  yNPpt1= NP1->GetBinContent(0);
  yNPpt2= NP2->GetBinContent(0);
  yNPpt3= NP3->GetBinContent(0);
  yNPpt4= NP4->GetBinContent(0);
  yNPpt5= NP5->GetBinContent(0);
  yNPpt1Err= NP1->GetBinError(0);
  yNPpt2Err= NP2->GetBinError(0);
  yNPpt3Err= NP3->GetBinError(0);
  yNPpt4Err= NP4->GetBinError(0);
  yNPpt5Err= NP5->GetBinError(0);
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

  gNPLow->SetPoint(0,(ptBin[0]+ptBin[1])/2,yNPpt1);
  gNPLow->SetPointError(0,0,yNPpt1Err);
  gNP->SetPoint(0,(ptBin[1]+ptBin[2])/2,yNPpt2);
  gNP->SetPointError(0,0,yNPpt2Err);
  gNP->SetPoint(1,(ptBin[2]+ptBin[3])/2,yNPpt3);
  gNP->SetPointError(1,0,yNPpt3Err);
  gNP->SetPoint(2,(ptBin[3]+ptBin[4])/2,yNPpt4);
  gNP->SetPointError(2,0,yNPpt4Err);
  gNP->SetPoint(3,(ptBin[4]+ptBin[5])/2,yNPpt5);
  gNP->SetPointError(3,0,yNPpt5Err);

  gPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gPR->GetXaxis()->CenterTitle();
  gPR->GetYaxis()->SetTitle("v_{3}");
  gPR->GetYaxis()->CenterTitle();
  gPR->GetXaxis()->SetLimits(0.,50);
  gPR->SetMinimum(-0.08);
  gPR->SetMaximum(0.20);

  gPR->SetLineColor(kMagenta+3);
  gPR->SetMarkerColor(kMagenta+3);
  gPR->SetMarkerSize(2);
  gPRLow->SetLineColor(kBlue+2);
  gPRLow->SetMarkerColor(kBlue+2);
  gPRLow->SetMarkerSize(2);

  gNP->SetLineColor(kOrange-5);
  gNP->SetMarkerColor(kOrange-5);
  gNP->SetMarkerSize(3);
  gNP->SetMarkerStyle(33);

  gNPLow->SetLineColor(kRed+2);
  gNPLow->SetMarkerColor(kRed+2);
  gNPLow->SetMarkerSize(3);
  gNPLow->SetMarkerStyle(33);


  float pos_x = 0.15;
  float pos_x_mass = 0.25;
  float pos_y = 0.90;
  float pos_y_diff = 0.071;
  int text_color = 1;
  float text_size = 25;
  TString perc = "%";

  TLegend *leg = new TLegend(0.58,0.6,0.8,0.8);
  SetLegendStyle(leg);
  leg->AddEntry(gPRLow,"Prompt J/#psi, 1.6 < |y^{#mu#mu}| < 2.4","pl");
  leg->AddEntry(gNPLow,"Nonprompt J/#psi,1.6 < |y^{#mu#mu}| < 2.4","pl");
  leg->AddEntry(gPR,"Prompt J/#psi, |y^{#mu#mu}| < 2.4","pl");
  leg->AddEntry(gNP,"Nonprompt J/#psi, |y^{#mu#mu}| < 2.4","pl");


  TCanvas *c1 = new TCanvas("c1","",900,800);
  c1->cd();
  gPR->Draw("APE");
  gNP->Draw("pe");
  gPRLow->Draw("pe");
  gNPLow->Draw("pe");
  leg->Draw("sames");
  jumSun(0,0,90,0);

  drawText("Scalar Product", pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  //drawText("6.5 - 50 GeV/c, |y^{#mu#mu}| < 2.4, SP", pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  drawText("Cent. 10-60 %", pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);

  //leg = new TLegend(0.74, 0.70, 0.925, 0.77);
  //SetLegendStyle(leg);
  //leg->SetTextSize(0.044);
  //leg->AddEntry(hist_v2,"J/#psi","pe");
  //leg->Draw("same");
  //globtex->DrawLatex(0.23, sz_init, "1.6 < |y| < 2.4");
  //globtex->DrawLatex(0.23, sz_init-sz_step, Form("Cent. 10-60%s",perc.Data()));
  CMS_lumi_v2mass(c1,iPeriod,iPos);
  c1->Update();
  c1->SaveAs(Form("new_v3_pT.pdf"));
  c1->SaveAs(Form("new_v3_pT.png"));
  //c1->SaveAs("v2_pt_Cent1060.pdf");
  //TFile* fSys2S = new TFile("../Systematic/merged_sys.root","read");
  //TH1D* hPt1S = (TH1D*) fSys->Get("hSys_pt_merged");
  //TH1D* hPt2S = (TH1D*) fSys2S->Get("hSys_pt_merged");
  //
  //double yPt1Ssys = hPt1S->GetBinContent(1);
  //double yPt2Ssys = hPt2S->GetBinContent(1);

  //cout << "Ptegr1ated 1S : " << Form("%.3f",yPt1S) << " +/- " << Form("%.3f",yPt1SErr) << " +/- " << Form("%.3f",yPt1Ssys) << endl;
  //cout << "Ptegr1ated 2S : " << Form("%.3f",yPt2S) << " +/- " << Form("%.3f",yPt2SErr) << " +/- " << Form("%.3f",yPt2Ssys) << endl;

  TFile *out = new TFile(Form("out_new_v3_vs_pT.root"),"recreate");
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
