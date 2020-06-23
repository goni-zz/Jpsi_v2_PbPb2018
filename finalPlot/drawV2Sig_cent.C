#include <iostream>
#include "../Style.h"
#include "../commonUtility.h"
#include "../cutsAndBin.h"
//#include "../CMS_lumi_v2mass.C"
//#include "../CMS_lumi_v2mass.h"
#include "../CMS_lumi_square.C"
#include "../tdrstyle.C"

using namespace std;

void drawV2Sig_cent(int PR=0){

  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
 
  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TGraphErrors* gr = new TGraphErrors();
  TString bCont;
  if(PR==0) bCont="Prompt";
  else if(PR==1) bCont="NonPrompt";
                                                                                             
  if(PR==0){
  TFile *fPRpt1 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality0-20_m2.6-3.5_OS.root","read");
  TFile *fPRpt2 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality20-40_m2.6-3.5_OS.root","read"); 
  TFile *fPRpt3 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality40-60_m2.6-3.5_OS.root","read"); 
  TFile *fPRpt4 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality60-80_m2.6-3.5_OS.root","read"); 
  TFile *fPRpt5 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality80-100_m2.6-3.5_OS.root","read");
  TFile *fPRpt6 = new TFile("../roots/v2mass_fit/SimFitResult_Prompt_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality100-180_m2.6-3.5_OS.root","read");

  TGraphErrors* gPRpt1 = (TGraphErrors*) fPRpt1 -> Get("v2vspt");  
  TGraphErrors* gPRpt2 = (TGraphErrors*) fPRpt2 -> Get("v2vspt");  
  TGraphErrors* gPRpt3 = (TGraphErrors*) fPRpt3 -> Get("v2vspt");  
  TGraphErrors* gPRpt4 = (TGraphErrors*) fPRpt4 -> Get("v2vspt");  
  TGraphErrors* gPRpt5 = (TGraphErrors*) fPRpt5 -> Get("v2vspt");  
  TGraphErrors* gPRpt6 = (TGraphErrors*) fPRpt6 -> Get("v2vspt");  

  double xPRpt1, yPRpt1, xPRpt2, yPRpt2, xPRpt3, yPRpt3, xPRpt4, yPRpt4, xPRpt5, yPRpt5, xPRpt6, yPRpt6;
  double xPRpt1Err, yPRpt1Err, xPRpt2Err, yPRpt2Err, xPRpt3Err, yPRpt3Err, xPRpt4Err, yPRpt4Err, xPRpt5Err, yPRpt5Err, xPRpt6Err, yPRpt6Err;
 
  gPRpt1->GetPoint(0, xPRpt1, yPRpt1);
  gPRpt2->GetPoint(0, xPRpt2, yPRpt2);
  gPRpt3->GetPoint(0, xPRpt3, yPRpt3);
  gPRpt4->GetPoint(0, xPRpt4, yPRpt4);
  gPRpt5->GetPoint(0, xPRpt5, yPRpt5);
  gPRpt6->GetPoint(0, xPRpt6, yPRpt6);
  yPRpt1Err= gPRpt1->GetErrorY(0);
  yPRpt2Err= gPRpt2->GetErrorY(0);
  yPRpt3Err= gPRpt3->GetErrorY(0);
  yPRpt4Err= gPRpt4->GetErrorY(0);
  yPRpt5Err= gPRpt5->GetErrorY(0);
  yPRpt6Err= gPRpt6->GetErrorY(0);
  
  gr->SetPoint(0,5,yPRpt1);
  gr->SetPointError(0,5,yPRpt1Err);
  gr->SetPoint(1,15,yPRpt2);
  gr->SetPointError(1,5,yPRpt2Err);
  gr->SetPoint(2,25,yPRpt3);
  gr->SetPointError(2,5,yPRpt3Err);
  gr->SetPoint(3,35,yPRpt4);
  gr->SetPointError(3,5,yPRpt4Err);
  gr->SetPoint(4,45,yPRpt5);
  gr->SetPointError(4,5,yPRpt5Err);
  gr->SetPoint(5,70,yPRpt6);
  gr->SetPointError(5,20,yPRpt6Err);

  }
  else if(PR==1){
  TFile *fNPpt1 = new TFile("../roots/v2mass_fit/SimFitResult_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality0-20_m2.6-3.5_OS.root","read");   
  TFile *fNPpt2 = new TFile("../roots/v2mass_fit/SimFitResult_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality20-60_m2.6-3.5_OS.root","read");  
  TFile *fNPpt3 = new TFile("../roots/v2mass_fit/SimFitResult_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality60-100_m2.6-3.5_OS.root","read"); 
  TFile *fNPpt4 = new TFile("../roots/v2mass_fit/SimFitResult_pt3.0-50.0_y1.6-2.4_muPt0.0_centrality100-180_m2.6-3.5_OS.root","read");

  TGraphErrors* gNPpt1 = (TGraphErrors*) fNPpt1 -> Get("v2vspt");  
  TGraphErrors* gNPpt2 = (TGraphErrors*) fNPpt2 -> Get("v2vspt");  
  TGraphErrors* gNPpt3 = (TGraphErrors*) fNPpt3 -> Get("v2vspt");  
  TGraphErrors* gNPpt4 = (TGraphErrors*) fNPpt4 -> Get("v2vspt");  
  
  double xNPpt1, yNPpt1, xNPpt2, yNPpt2, xNPpt3, yNPpt3, xNPpt4, yNPpt4;
  double xNPpt1Err, yNPpt1Err, xNPpt2Err, yNPpt2Err, xNPpt3Err, yNPpt3Err, xNPpt4Err, yNPpt4Err;

  gNPpt1->GetPoint(0, xNPpt1, yNPpt1);
  gNPpt2->GetPoint(0, xNPpt2, yNPpt2);
  gNPpt3->GetPoint(0, xNPpt3, yNPpt3);
  gNPpt4->GetPoint(0, xNPpt4, yNPpt4);
  yNPpt1Err= gNPpt1->GetErrorY(0);
  yNPpt2Err= gNPpt2->GetErrorY(0);
  yNPpt3Err= gNPpt3->GetErrorY(0);
  yNPpt4Err= gNPpt4->GetErrorY(0);

  gr->SetPoint(0,5+0.1,yNPpt1);
  gr->SetPointError(0,5,yNPpt1Err);
  gr->SetPoint(1,double(10+30)/2,yNPpt2);
  gr->SetPointError(1,double(30-10)/2,yNPpt2Err);
  gr->SetPoint(2,double(30+50)/2,yNPpt3);
  gr->SetPointError(2,double(50-30)/2,yNPpt3Err);
  gr->SetPoint(3,double(50+90)/2+0.1,yNPpt4);
  gr->SetPointError(3,double(90-50)/2,yNPpt4Err);
  }

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

  if(PR==0){
  gr->SetLineColor(kMagenta+3);
  gr->SetMarkerColor(kMagenta+3);
  gr->SetMarkerSize(2);}
  else if(PR==1){
  gr->SetLineColor(kOrange-5);
  gr->SetMarkerColor(kOrange-5);
  gr->SetMarkerSize(3);
  gr->SetMarkerStyle(33);}

  float pos_x = 0.23;
  float pos_x_mass = 0.25;
  float pos_y = 0.90;
  float pos_y_diff = 0.071;
  int text_color = 1;
  float text_size = 25;
  TString perc = "%";

  TLegend *leg = new TLegend(0.7,0.65,0.9,0.8);
  SetLegendStyle(leg);
  if(PR==0){
  leg->AddEntry(gr,"Prompt J/#psi","p");}
  else if(PR==1){
  leg->AddEntry(gr,"Nonprompt J/#psi","p");}

  TCanvas *c1 = new TCanvas("c1","",900,800);
  c1->cd();
  gr->Draw("AP");
  leg->Draw("sames");
  jumSun(0,0,100,0);

  drawText("3 < p_{T} < 50 GeV/c, SP", pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  drawText("1.6 < |y^{#mu#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);

  //leg = new TLegend(0.74, 0.70, 0.925, 0.77);
  //SetLegendStyle(leg);
  //leg->SetTextSize(0.044);
  //leg->AddEntry(hist_v2,"J/#psi","pe");
  //leg->Draw("same");
  //globtex->DrawLatex(0.23, sz_init, "1.6 < |y| < 2.4");
  //globtex->DrawLatex(0.23, sz_init-sz_step, Form("Cent. 10-60%s",perc.Data()));
  CMS_lumi_square(c1,iPeriod,iPos);
  c1->Update();
  c1->SaveAs(Form("%s_v2_cent.pdf",bCont.Data()));
  //c1->SaveAs("v2_pt_Cent1060.pdf");
  //TFile* fSys2S = new TFile("../Systematic/merged_sys.root","read");
  //TH1D* hPt1S = (TH1D*) fSys->Get("hSys_pt_merged");
  //TH1D* hPt2S = (TH1D*) fSys2S->Get("hSys_pt_merged");
  //
  //double yPt1Ssys = hPt1S->GetBinContent(1);
  //double yPt2Ssys = hPt2S->GetBinContent(1);

  //cout << "Ptegr1ated 1S : " << Form("%.3f",yPt1S) << " +/- " << Form("%.3f",yPt1SErr) << " +/- " << Form("%.3f",yPt1Ssys) << endl;
  //cout << "Ptegr1ated 2S : " << Form("%.3f",yPt2S) << " +/- " << Form("%.3f",yPt2SErr) << " +/- " << Form("%.3f",yPt2Ssys) << endl;

  TFile *out = new TFile(Form("out_%s_v2_vs_cent.root",bCont.Data()),"recreate");
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
