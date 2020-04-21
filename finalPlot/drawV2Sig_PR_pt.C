#include <iostream>
#include "../Style.h"
#include "../commonUtility.h"
#include "../cutsAndBin.h"
//#include "../CMS_lumi_v2mass.C"
//#include "../CMS_lumi_v2mass.h"
#include "../CMS_lumi_square.C"
#include "../tdrstyle.C"

using namespace std;

void drawV2Sig_PR_pt(){

  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
 
  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile *fPRpt1 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-4.0_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS.root","read");
  TFile *fPRpt2 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt4.0-5.0_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS.root","read");
  TFile *fPRpt3 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt5.0-6.5_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS.root","read");
  TFile *fPRpt4 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt6.5-8.0_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS.root","read");
  TFile *fPRpt5 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt8.0-12.0_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS.root","read");
  TFile *fPRpt6 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt12.0-30.0_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS.root","read");

  TFile *fNPpt1 = new TFile("../roots/v2mass_fit/SimFitResult_pt3.0-4.5_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS.root","read");
  TFile *fNPpt2 = new TFile("../roots/v2mass_fit/SimFitResult_pt4.5-6.5_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS.root","read");
  TFile *fNPpt3 = new TFile("../roots/v2mass_fit/SimFitResult_pt6.5-9.0_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS.root","read");
  TFile *fNPpt4 = new TFile("../roots/v2mass_fit/SimFitResult_pt9.0-12.0_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS.root","read");
  TFile *fNPpt5 = new TFile("../roots/v2mass_fit/SimFitResult_pt12.0-30.0_y1.6-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS.root","read"); 
  
  TGraphErrors* gPRpt1 = (TGraphErrors*) fPRpt1 -> Get("v2vspt");  
  TGraphErrors* gPRpt2 = (TGraphErrors*) fPRpt2 -> Get("v2vspt");  
  TGraphErrors* gPRpt3 = (TGraphErrors*) fPRpt3 -> Get("v2vspt");  
  TGraphErrors* gPRpt4 = (TGraphErrors*) fPRpt4 -> Get("v2vspt");  
  TGraphErrors* gPRpt5 = (TGraphErrors*) fPRpt5 -> Get("v2vspt");  
  TGraphErrors* gPRpt6 = (TGraphErrors*) fPRpt6 -> Get("v2vspt");  

  TGraphErrors* gNPpt1 = (TGraphErrors*) fNPpt1 -> Get("v2vspt");  
  TGraphErrors* gNPpt2 = (TGraphErrors*) fNPpt2 -> Get("v2vspt");  
  TGraphErrors* gNPpt3 = (TGraphErrors*) fNPpt3 -> Get("v2vspt");  
  TGraphErrors* gNPpt4 = (TGraphErrors*) fNPpt4 -> Get("v2vspt");  
  TGraphErrors* gNPpt5 = (TGraphErrors*) fNPpt5 -> Get("v2vspt");  

  double xPRpt1, yPRpt1, xPRpt2, yPRpt2, xPRpt3, yPRpt3, xPRpt4, yPRpt4, xPRpt5, yPRpt5, xPRpt6, yPRpt6;
  double xPRpt1Err, yPRpt1Err, xPRpt2Err, yPRpt2Err, xPRpt3Err, yPRpt3Err, xPRpt4Err, yPRpt4Err, xPRpt5Err, yPRpt5Err, xPRpt6Err, yPRpt6Err;
  
  double xNPpt1, yNPpt1, xNPpt2, yNPpt2, xNPpt3, yNPpt3, xNPpt4, yNPpt4, xNPpt5, yNPpt5;
  double xNPpt1Err, yNPpt1Err, xNPpt2Err, yNPpt2Err, xNPpt3Err, yNPpt3Err, xNPpt4Err, yNPpt4Err, xNPpt5Err, yNPpt5Err;

  gPRpt1->GetPoint(0, xPRpt1, yPRpt1);
  gPRpt2->GetPoint(0, xPRpt2, yPRpt2);
  gPRpt3->GetPoint(0, xPRpt3, yPRpt3);
  gPRpt4->GetPoint(0, xPRpt4, yPRpt4);
  gPRpt5->GetPoint(0, xPRpt5, yPRpt5);
  gPRpt6->GetPoint(0, xPRpt6, yPRpt6);
  gNPpt1->GetPoint(0, xNPpt1, yNPpt1);
  gNPpt2->GetPoint(0, xNPpt2, yNPpt2);
  gNPpt3->GetPoint(0, xNPpt3, yNPpt3);
  gNPpt4->GetPoint(0, xNPpt4, yNPpt4);
  gNPpt5->GetPoint(0, xNPpt5, yNPpt5);
  yPRpt1Err= gPRpt1->GetErrorY(0);
  yPRpt2Err= gPRpt2->GetErrorY(0);
  yPRpt3Err= gPRpt3->GetErrorY(0);
  yPRpt4Err= gPRpt4->GetErrorY(0);
  yPRpt5Err= gPRpt5->GetErrorY(0);
  yPRpt6Err= gPRpt6->GetErrorY(0);
  yNPpt1Err= gNPpt1->GetErrorY(0);
  yNPpt2Err= gNPpt2->GetErrorY(0);
  yNPpt3Err= gNPpt3->GetErrorY(0);
  yNPpt4Err= gNPpt4->GetErrorY(0);
  yNPpt5Err= gNPpt5->GetErrorY(0);
  TGraphErrors* gr1 = new TGraphErrors();
  TGraphErrors* gr2 = new TGraphErrors();

  gr1->SetPoint(0,3.5,yPRpt1);
  gr1->SetPointError(0,0.5,yPRpt1Err);
  gr1->SetPoint(1,4.5,yPRpt2);
  gr1->SetPointError(1,0.5,yPRpt2Err);
  gr1->SetPoint(2,(5+6.5)/2,yPRpt3);
  gr1->SetPointError(2,(6.5-5)/2,yPRpt3Err);
  gr1->SetPoint(3,(6.5+8)/2,yPRpt4);
  gr1->SetPointError(3,(8-6.5)/2,yPRpt4Err);
  gr1->SetPoint(4,(8+12)/2,yPRpt5);
  gr1->SetPointError(4,(12-8)/2,yPRpt5Err);
  gr1->SetPoint(5,(12+30)/2,yPRpt6);
  gr1->SetPointError(5,(30-12)/2,yPRpt6Err);

  gr2->SetPoint(0,(3+4.5)/2,yNPpt1);
  gr2->SetPointError(0,(4.5-3)/2,yNPpt1Err);
  gr2->SetPoint(1,(4.5+6.5)/2,yNPpt2);
  gr2->SetPointError(1,(6.5-4.5)/2,yNPpt2Err);
  gr2->SetPoint(2,(6.5+9)/2,yNPpt3);
  gr2->SetPointError(2,(9-6.5)/2,yNPpt3Err);
  gr2->SetPoint(3,double(9+12)/2,yNPpt4);
  gr2->SetPointError(3,double(12-9)/2,yNPpt4Err);
  gr2->SetPoint(4,(12+30)/2,yNPpt5);
  gr2->SetPointError(4,(30-12)/2,yNPpt5Err);
  TMultiGraph *mg = new TMultiGraph();
  //mg->SetTitle("Exclusion graphs");
  mg->Add(gr1);
  mg->Add(gr2);
  mg->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetTitle("v_{2}(s)");
  mg->GetYaxis()->CenterTitle();
  mg->GetXaxis()->SetLimits(0.,30);
  mg->SetMinimum(-0.08);
  mg->SetMaximum(0.25);

  //// for cent
  mg->GetXaxis()->SetTitleSize(0.06*1.0);
  mg->GetYaxis()->SetTitleSize(0.06*1.0);
  mg->GetXaxis()->SetLabelSize(0.05*1.0);
  mg->GetYaxis()->SetLabelSize(0.05*1.0);

  gr1->SetLineColor(kMagenta+3);
  gr1->SetMarkerColor(kMagenta+3);
  gr1->SetMarkerSize(1.7);
  gr2->SetLineColor(kOrange-5);
  gr2->SetMarkerColor(kOrange-5);
  gr2->SetMarkerSize(2.2);
  gr2->SetMarkerStyle(33);

  float pos_x = 0.23;
  float pos_x_mass = 0.25;
  float pos_y = 0.90;
  float pos_y_diff = 0.071;
  int text_color = 1;
  float text_size = 25;
  TString perc = "%";

  TLegend *leg = new TLegend(0.7,0.65,0.9,0.8);
  SetLegendStyle(leg);
  leg->AddEntry(gr1,"Prompt J/#psi","p");
  leg->AddEntry(gr2,"Nonprompt J/#psi","p");

  TCanvas *c1 = new TCanvas("c1","",900,800);
  c1->cd();
  mg->Draw("AP");
  //gr1->Draw("AP");
  //gr2->Draw("AP");
  leg->Draw("sames");
  jumSun(0,0,30,0);

  drawText("1.6 < |y^{#mu#mu}| < 2.4, SP", pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  drawText("Cent.10-60%", pos_x_mass,pos_y-pos_y_diff*1.5,text_color,text_size);

  //leg = new TLegend(0.74, 0.70, 0.925, 0.77);
  //SetLegendStyle(leg);
  //leg->SetTextSize(0.044);
  //leg->AddEntry(hist_v2,"J/#psi","pe");
  //leg->Draw("same");
  //globtex->DrawLatex(0.23, sz_init, "1.6 < |y| < 2.4");
  //globtex->DrawLatex(0.23, sz_init-sz_step, Form("Cent. 10-60%s",perc.Data()));
  CMS_lumi_square(c1,iPeriod,iPos);
  c1->Update();
  c1->SaveAs("Prompt_v2_pt.pdf");
  //c1->SaveAs("v2_pt_Cent1060.pdf");
  //TFile* fSys2S = new TFile("../Systematic/merged_sys.root","read");
  //TH1D* hPt1S = (TH1D*) fSys->Get("hSys_pt_merged");
  //TH1D* hPt2S = (TH1D*) fSys2S->Get("hSys_pt_merged");
  //
  //double yPt1Ssys = hPt1S->GetBinContent(1);
  //double yPt2Ssys = hPt2S->GetBinContent(1);

  //cout << "Ptegr1ated 1S : " << Form("%.3f",yPt1S) << " +/- " << Form("%.3f",yPt1SErr) << " +/- " << Form("%.3f",yPt1Ssys) << endl;
  //cout << "Ptegr1ated 2S : " << Form("%.3f",yPt2S) << " +/- " << Form("%.3f",yPt2SErr) << " +/- " << Form("%.3f",yPt2Ssys) << endl;

  TFile *out = new TFile("out_Prompt_v2_vs_pt.root","recreate");
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
