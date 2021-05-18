#include <iostream>
#include "../Style.h"
#include "../commonUtility.h"
#include "../cutsAndBin.h"
#include "../CMS_lumi_v2mass.C"
//#include "../CMS_lumi_v2mass.h"
//#include "../CMS_lumi_square.C"
#include "../tdrstyle.C"

using namespace std;

void drawV2Sig_pt_new(){

  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;

  const int nPtBins=6;
  float ptBin[nPtBins+1]={3,4.5,6.5,7.5,9,12,50};
 
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

  TFile *fPRpt1 = new TFile("v2_pt3.0-4.5_y1.6-2.4_muPt0.0_centrality40-80.root","read");
  TFile *fPRpt2 = new TFile("v2_pt4.5-6.5_y1.6-2.4_muPt0.0_centrality40-80.root","read");
  TFile *fPRpt3 = new TFile("v2_pt6.5-7.5_y0.0-2.4_muPt0.0_centrality40-80.root","read");
  TFile *fPRpt4 = new TFile("v2_pt7.5-9.0_y0.0-2.4_muPt0.0_centrality40-80.root","read");
  TFile *fPRpt5 = new TFile("v2_pt9.0-12.0_y0.0-2.4_muPt0.0_centrality40-80.root","read");
  TFile *fPRpt6 = new TFile("v2_pt12.0-50.0_y0.0-2.4_muPt0.0_centrality40-80.root","read");

  TH1D * PRpt1 = (TH1D*) fPRpt1 -> Get("prv2"); 
  TH1D * PRpt2 = (TH1D*) fPRpt2 -> Get("prv2"); 
  TH1D * PRpt3 = (TH1D*) fPRpt3 -> Get("prv2"); 
  TH1D * PRpt4 = (TH1D*) fPRpt4 -> Get("prv2"); 
  TH1D * PRpt5 = (TH1D*) fPRpt5 -> Get("prv2"); 
  TH1D * PRpt6 = (TH1D*) fPRpt6 -> Get("prv2"); 

  TH1D * NPpt1 = (TH1D*) fPRpt1 -> Get("npv2"); 
  TH1D * NPpt2 = (TH1D*) fPRpt2 -> Get("npv2"); 
  TH1D * NPpt3 = (TH1D*) fPRpt3 -> Get("npv2"); 
  TH1D * NPpt4 = (TH1D*) fPRpt4 -> Get("npv2"); 
  TH1D * NPpt5 = (TH1D*) fPRpt5 -> Get("npv2"); 
  TH1D * NPpt6 = (TH1D*) fPRpt6 -> Get("npv2"); 

  double xPRpt1, yPRpt1, xPRpt2, yPRpt2, xPRpt3, yPRpt3, xPRpt4, yPRpt4, xPRpt5, yPRpt5, xPRpt6, yPRpt6;
  double xPRpt1Err, yPRpt1Err, xPRpt2Err, yPRpt2Err, xPRpt3Err, yPRpt3Err, xPRpt4Err, yPRpt4Err, xPRpt5Err, yPRpt5Err, xPRpt6Err, yPRpt6Err;
 
  yPRpt1= PRpt1->GetBinContent(0);
  yPRpt2= PRpt2->GetBinContent(0);
  yPRpt3= PRpt3->GetBinContent(0);
  yPRpt4= PRpt4->GetBinContent(0);
  yPRpt5= PRpt5->GetBinContent(0);
  yPRpt6= PRpt6->GetBinContent(0);
  yPRpt1Err= PRpt1->GetBinError(0);
  yPRpt2Err= PRpt2->GetBinError(0);
  yPRpt3Err= PRpt3->GetBinError(0);
  yPRpt4Err= PRpt4->GetBinError(0);
  yPRpt5Err= PRpt5->GetBinError(0);
  yPRpt6Err= PRpt6->GetBinError(0);
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
  gPRLow->SetPoint(1,(ptBin[1]+ptBin[2])/2,yPRpt2);
  gPRLow->SetPointError(1,0,yPRpt2Err);
  gPR->SetPoint(0,(ptBin[2]+ptBin[3])/2,yPRpt3);
  gPR->SetPointError(0,0,yPRpt3Err);
  gPR->SetPoint(1,(ptBin[3]+ptBin[4])/2,yPRpt4);
  gPR->SetPointError(1,0,yPRpt4Err);
  gPR->SetPoint(2,(ptBin[4]+ptBin[5])/2,yPRpt5);
  gPR->SetPointError(2,0,yPRpt5Err);
  gPR->SetPoint(3,(ptBin[5]+ptBin[6])/2,yPRpt6);
  gPR->SetPointError(3,0,yPRpt6Err);

  double xNPpt1, yNPpt1, xNPpt2, yNPpt2, xNPpt3, yNPpt3, xNPpt4, yNPpt4, xNPpt5, yNPpt5, xNPpt6, yNPpt6;
  double xNPpt1Err, yNPpt1Err, xNPpt2Err, yNPpt2Err, xNPpt3Err, yNPpt3Err, xNPpt4Err, yNPpt4Err, xNPpt5Err, yNPpt5Err, xNPpt6Err, yNPpt6Err;
 
  yNPpt1= NPpt1->GetBinContent(0);
  yNPpt2= NPpt2->GetBinContent(0);
  yNPpt3= NPpt3->GetBinContent(0);
  yNPpt4= NPpt4->GetBinContent(0);
  yNPpt5= NPpt5->GetBinContent(0);
  yNPpt6= NPpt6->GetBinContent(0);
  yNPpt1Err= NPpt1->GetBinError(0);
  yNPpt2Err= NPpt2->GetBinError(0);
  yNPpt3Err= NPpt3->GetBinError(0);
  yNPpt4Err= NPpt4->GetBinError(0);
  yNPpt5Err= NPpt5->GetBinError(0);
  yNPpt6Err= NPpt6->GetBinError(0);
//  
  //gNPLow->SetPoint(0,(ptBin[0]+ptBin[1])/2,yNPpt1);
  //gNPLow->SetPointError(0,(ptBin[1]-ptBin[0])/2,yNPpt1Err);
  //gNPLow->SetPoint(1,(ptBin[1]+ptBin[2])/2,yNPpt2);
  //gNPLow->SetPointError(1,(ptBin[2]-ptBin[1])/2,yNPpt2Err);
  //gNP->SetPoint(0,(ptBin[2]+ptBin[3])/2,yNPpt3);
  //gNP->SetPointError(0,(ptBin[3]-ptBin[2])/2,yNPpt3Err);
  //gNP->SetPoint(1,(ptBin[3]+ptBin[4])/2,yNPpt4);
  //gNP->SetPointError(1,(ptBin[4]-ptBin[3])/2,yNPpt4Err);
  //gNP->SetPoint(2,(ptBin[4]+ptBin[5])/2,yNPpt5);
  //gNP->SetPointError(2,(ptBin[5]-ptBin[4])/2,yNPpt5Err);
  //gNP->SetPoint(3,(ptBin[5]+ptBin[6])/2,yNPpt6);
  //gNP->SetPointError(3,(ptBin[6]-ptBin[5])/2,yNPpt6Err);
  gNPLow->SetPoint(0,(ptBin[0]+ptBin[1])/2,yNPpt1);
  gNPLow->SetPointError(0,0,yNPpt1Err);
  gNPLow->SetPoint(1,(ptBin[1]+ptBin[2])/2,yNPpt2);
  gNPLow->SetPointError(1,0,yNPpt2Err);
  gNP->SetPoint(0,(ptBin[2]+ptBin[3])/2,yNPpt3);
  gNP->SetPointError(0,0,yNPpt3Err);
  gNP->SetPoint(1,(ptBin[3]+ptBin[4])/2,yNPpt4);
  gNP->SetPointError(1,0,yNPpt4Err);
  gNP->SetPoint(2,(ptBin[4]+ptBin[5])/2,yNPpt5);
  gNP->SetPointError(2,0,yNPpt5Err);
  gNP->SetPoint(3,(ptBin[5]+ptBin[6])/2,yNPpt6);
  gNP->SetPointError(3,0,yNPpt6Err);

  gPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gPR->GetXaxis()->CenterTitle();
  gPR->GetYaxis()->SetTitle("v_{2}");
  gPR->GetYaxis()->CenterTitle();
  gPR->GetXaxis()->SetLimits(0.,50);
  gPR->SetMinimum(-0.08);
  gPR->SetMaximum(0.25);

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

  TLegend *leg = new TLegend(0.58,0.55,0.8,0.8);
  SetLegendStyle(leg);
  leg->AddEntry(gPRLow,"Prompt J/#psi, 1.6 < |y^{#mu#mu}| < 2.4","p");
  leg->AddEntry(gNPLow,"Nonprompt J/#psi,1.6 < |y^{#mu#mu}| < 2.4","p");
  leg->AddEntry(gPR,"Prompt J/#psi, |y^{#mu#mu}| < 2.4","p");
  leg->AddEntry(gNP,"Nonprompt J/#psi, |y^{#mu#mu}| < 2.4","p");

  TCanvas *c1 = new TCanvas("c1","",900,800);
  c1->cd();
  gPR->Draw("AP");
  gNP->Draw("p");
  gPRLow->Draw("p");
  gNPLow->Draw("p");
  leg->Draw("sames");
  jumSun(0,0,50,0);

  drawText("3 - 50 GeV/c, SP", pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  //drawText("6.5 - 50 GeV/c, |y^{#mu#mu}| < 2.4, SP", pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  drawText("Cent.20-40%", pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);

  //leg = new TLegend(0.74, 0.70, 0.925, 0.77);
  //SetLegendStyle(leg);
  //leg->SetTextSize(0.044);
  //leg->AddEntry(hist_v2,"J/#psi","pe");
  //leg->Draw("same");
  //globtex->DrawLatex(0.23, sz_init, "1.6 < |y| < 2.4");
  //globtex->DrawLatex(0.23, sz_init-sz_step, Form("Cent. 10-60%s",perc.Data()));
  CMS_lumi_v2mass(c1,iPeriod,iPos);
  c1->Update();
  c1->SaveAs(Form("new_v2_pt_cent_20_40.pdf"));
  //c1->SaveAs("v2_pt_Cent1060.pdf");
  //TFile* fSys2S = new TFile("../Systematic/merged_sys.root","read");
  //TH1D* hPt1S = (TH1D*) fSys->Get("hSys_pt_merged");
  //TH1D* hPt2S = (TH1D*) fSys2S->Get("hSys_pt_merged");
  //
  //double yPt1Ssys = hPt1S->GetBinContent(1);
  //double yPt2Ssys = hPt2S->GetBinContent(1);

  //cout << "Ptegr1ated 1S : " << Form("%.3f",yPt1S) << " +/- " << Form("%.3f",yPt1SErr) << " +/- " << Form("%.3f",yPt1Ssys) << endl;
  //cout << "Ptegr1ated 2S : " << Form("%.3f",yPt2S) << " +/- " << Form("%.3f",yPt2SErr) << " +/- " << Form("%.3f",yPt2Ssys) << endl;

  TFile *out = new TFile(Form("out_new_v2_vs_pt_20_40.root"),"recreate");
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
