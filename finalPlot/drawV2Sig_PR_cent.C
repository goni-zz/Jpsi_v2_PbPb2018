#include <iostream>
#include "../Style.h"
#include "../commonUtility.h"
#include "../cutsAndBin.h"
#include "../CMS_lumi_v2mass.C"
#include "../tdrstyle.C"

using namespace std;

void drawV2Sig_PR_cent(){

  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
 
  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile *fPt1 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality0-20_m2.6-3.5_OS.root","read");
  TFile *fPt2 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality20-40_m2.6-3.5_OS.root","read");
  TFile *fPt3 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality40-60_m2.6-3.5_OS.root","read");
  TFile *fPt4 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality60-80_m2.6-3.5_OS.root","read");
  TFile *fPt5 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality80-100_m2.6-3.5_OS.root","read");
  TFile *fPt6 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality100-180_m2.6-3.5_OS.root","read");
  
  TGraphErrors* gPt1 = (TGraphErrors*) fPt1 -> Get("v2vspt");  
  TGraphErrors* gPt2 = (TGraphErrors*) fPt2 -> Get("v2vspt");  
  TGraphErrors* gPt3 = (TGraphErrors*) fPt3 -> Get("v2vspt");  
  TGraphErrors* gPt4 = (TGraphErrors*) fPt4 -> Get("v2vspt");  
  TGraphErrors* gPt5 = (TGraphErrors*) fPt5 -> Get("v2vspt");  
  TGraphErrors* gPt6 = (TGraphErrors*) fPt6 -> Get("v2vspt");  
  double xPt1, yPt1, xPt2, yPt2, xPt3, yPt3, xPt4, yPt4, xPt5, yPt5, xPt6, yPt6;
  double xPt1Err, yPt1Err, xPt2Err, yPt2Err, xPt3Err, yPt3Err, xPt4Err, yPt4Err, xPt5Err, yPt5Err, xPt6Err, yPt6Err;
  
  gPt1->GetPoint(0, xPt1, yPt1);
  gPt2->GetPoint(0, xPt2, yPt2);
  gPt3->GetPoint(0, xPt3, yPt3);
  gPt4->GetPoint(0, xPt4, yPt4);
  gPt5->GetPoint(0, xPt5, yPt5);
  gPt6->GetPoint(0, xPt6, yPt6);

  yPt1Err= gPt1->GetErrorY(0);
  yPt2Err= gPt2->GetErrorY(0);
  yPt3Err= gPt3->GetErrorY(0);
  yPt4Err= gPt4->GetErrorY(0);
  yPt5Err= gPt5->GetErrorY(0);
  yPt6Err= gPt6->GetErrorY(0);

  TGraphErrors* gr = new TGraphErrors();

  gr->SetPoint(0,5,yPt1);
  gr->SetPointError(0,5,yPt1Err);
  gr->SetPoint(1,15,yPt2);
  gr->SetPointError(1,5,yPt2Err);
  gr->SetPoint(2,25,yPt3);
  gr->SetPointError(2,5,yPt3Err);
  gr->SetPoint(3,35,yPt4);
  gr->SetPointError(3,5,yPt4Err);
  gr->SetPoint(4,45,yPt5);
  gr->SetPointError(4,5,yPt5Err);
  gr->SetPoint(5,70,yPt6);
  gr->SetPointError(5,20,yPt6Err);

  gr->GetXaxis()->SetTitle("Centrality (%)");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle("v_{2}(s)");
  gr->GetYaxis()->CenterTitle();
  gr->GetXaxis()->SetLimits(0.,100);
  gr->SetMinimum(-0.08);
  gr->SetMaximum(0.25);

  //// for cent
  gr->GetXaxis()->SetTitleSize(0.06*1.0);
  gr->GetYaxis()->SetTitleSize(0.06*1.0);
  gr->GetXaxis()->SetLabelSize(0.05*1.0);
  gr->GetYaxis()->SetLabelSize(0.05*1.0);

  gr->SetLineColor(kMagenta+3);
  gr->SetMarkerColor(kMagenta+3);
  gr->SetMarkerSize(1.5);

  float pos_x = 0.23;
  float pos_x_mass = 0.25;
  float pos_y = 0.85;
  float pos_y_diff = 0.071;
  int text_color = 1;
  float text_size = 25;
  TString perc = "%";

  TLegend *leg = new TLegend(0.65,0.78,0.8,0.9);
  SetLegendStyle(leg);
  leg->AddEntry(gr,"Prompt J/#psi","p");

  TCanvas *c1 = new TCanvas("c1","",900,800);
  c1->cd();
  gr->Draw("AP");
  leg->Draw("same");
  jumSun(0,0,100,0);

  drawText("3 < p_{T} < 50 GeV/c", pos_x_mass,pos_y-pos_y_diff*0.5,text_color,text_size);
  drawText("1.6 < |y^{#mu#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff,text_color,text_size);

  c1->SaveAs("Prompt_v2_cent.pdf");

  //TFile* fSys2S = new TFile("../Systematic/merged_sys.root","read");
  //TH1D* hPt1S = (TH1D*) fSys->Get("hSys_pt_merged");
  //TH1D* hPt2S = (TH1D*) fSys2S->Get("hSys_pt_merged");
  //
  //double yPt1Ssys = hPt1S->GetBinContent(1);
  //double yPt2Ssys = hPt2S->GetBinContent(1);

  //cout << "Ptegrated 1S : " << Form("%.3f",yPt1S) << " +/- " << Form("%.3f",yPt1SErr) << " +/- " << Form("%.3f",yPt1Ssys) << endl;
  //cout << "Ptegrated 2S : " << Form("%.3f",yPt2S) << " +/- " << Form("%.3f",yPt2SErr) << " +/- " << Form("%.3f",yPt2Ssys) << endl;

  TFile *out = new TFile("out_Prompt_v2_vs_cent.root","recreate");
//  out->cd();
//  final_v2_int_g->SetName("gr_point_v2_vs_pt_Cent1060");
//  gsys_int->SetName("gr_sys_v2_vs_pt_Cent1060");
//  final_v2_int_g->Write();
//  gsys_int->Write(); 
//  for(int i=0;i<nCentBin;i++){
//    gsys[i]->Write();
//    final_v2_g[i]->Write();
//  }
//	return;
} 
