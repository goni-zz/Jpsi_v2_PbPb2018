#include <iostream>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include <TMath.h>
#include "../commonUtility.h"
#include "../cutsAndBin.h"
#include "../Style.h"
#include "../tdrstyle.C"
#include "../CMS_lumi_v2mass.C"
#include "../rootFitHeaders.h"

using namespace std;

void GetHistSqrt(TH1D* h1 =0, TH1D* h2=0);

void psuedo_proper_decay_length(int cLow = 120, int cHigh = 180,
    float ptLow =  3, float ptHigh = 30, 
    float yLow = 1.6, float yHigh = 2.4,
    float SiMuPtCut = 0, float massLow = 2.6, float massHigh =3.5, int PRMC=0, bool dimusign=true, bool fAccW = false, bool fEffW = false)
{

  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh) ;
  kineLabel = kineLabel + Form("_m%.1f-%.1f",massLow,massHigh);

  //Basic Setting
  gStyle->SetOptStat(0);
  //  setTDRStyle();

  TFile *fData;
  TFile *fPRMC;
  TFile *fNPMC;
  TTree *treeData;
  TTree *treePRMC;
  TTree *treeNPMC;

  fData  = new TFile("../skimmedFiles/OniaFlowSkim_JpsiTrig_AllPD_isMC0_HFNom_AddEP_200217.root");
  fPRMC  = new TFile("../skimmedFiles/OniaFlowSkim_Jpsi_MC_Prompt_RECO.root");
  fNPMC  = new TFile("../skimmedFiles/OniaFlowSkim_Jpsi_MC_NonPrompt_RECO.root");

  TCut mCut;
  TCut ptCut;
  TCut yCut;
  TCut cCut;
  TCut sglMuCut;
  TCut accCut;
  TCut dimuSign;
  TCut quality;
  quality = ("fabs(vz)<15");
  dimuSign = ("recoQQsign==0");
  mCut = ("mass>=2.6&&mass<=3.5");
  ptCut = Form("pt>%1.f&&pt<%.1f", ptLow, ptHigh);
  yCut  = Form("abs(y)>%.1f&&abs(y)<%.1f", yLow, yHigh);
  cCut  = Form("cBin>=%d&&cBin<=%d",cLow,cHigh);
  sglMuCut = Form("abs(eta1)>%.1f&&abs(eta1)<%.1f&&abs(eta2)>%.1f&&abs(eta2)<%.1f", yLow, yHigh, yLow, yHigh);
  accCut= ("( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || \
      ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || \
      ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )");

  cout<<kineLabel<<endl;
  TCut totalCut ="";
  totalCut = quality&&dimuSign&&mCut&&ptCut&&yCut&&cCut&&sglMuCut&&accCut;
  //totalCut = quality&&dimuSign&&mCut&&ptCut&&yCut&&cCut&&sglMuCut;
  cout<<totalCut<<endl;
  //if(forPR==0)consider="PR_eff";
  //else if(forPR==1)consider="NP_eff";

  float nbins=4000;
  float xmin =  -1;
  float xmax =   3;
  TH1D* h_mass = new TH1D("h_mass",";m_{#mu^{+}#mu^{-}};Counts/(0.02 GeV)",120,massLow,massHigh);
  TH1D* h_massCut = new TH1D("h_massCut",";m_{#mu^{+}#mu^{-}};Counts/(0.02 GeV)",120,massLow,massHigh);
  TH1D* h_decayData = new TH1D("h_decayData",";l_{J/#psi};Counts",nbins,xmin,xmax);
  TH1D* h_decayPRMC = new TH1D("h_decayPRMC",";l_{J/#psi};Counts",nbins,xmin,xmax);
  TH1D* h_decayNPMC = new TH1D("h_decayNPMC",";l_{J/#psi};Counts",nbins,xmin,xmax);
  TH1D* h_deffData = new TH1D("h_deffData",";l_{J/#psi};Efficiency",nbins,xmin,xmax);
  TH1D* h_deffPRMC = new TH1D("h_deffPRMC",";l_{J/#psi};Efficiency",nbins,xmin,xmax);
  TH1D* h_deffNPMC = new TH1D("h_deffNPMC",";l_{J/#psi};Efficiency",nbins,xmin,xmax);

  treeData = (TTree*)fData->Get("mmepevt");
  treePRMC = (TTree*)fPRMC->Get("mmepevt");
  treeNPMC = (TTree*)fNPMC->Get("mmepevt");

  treeData->Draw("ctau3D>>h_decayData",totalCut,"l");
  treePRMC->Draw("ctau3D>>h_decayPRMC",totalCut,"l");
  treeNPMC->Draw("ctau3D>>h_decayNPMC",totalCut,"l");

  cout<<"Data Tree: "<<treeData->GetEntries()<<", by total cut: "<<h_decayData->GetEntries()<<", "<<h_decayData->GetEntries()/treeData->GetEntries()*100<<"%"<<endl;
  cout<<"PRMC Tree: "<<treePRMC->GetEntries()<<", by total cut: "<<h_decayPRMC->GetEntries()<<", "<<h_decayData->GetEntries()/treePRMC->GetEntries()*100<<"%"<<endl;
  cout<<"NPMC Tree: "<<treeNPMC->GetEntries()<<", by total cut: "<<h_decayNPMC->GetEntries()<<", "<<h_decayData->GetEntries()/treeNPMC->GetEntries()*100<<"%"<<endl;

  TLine *lcutv;
  TLine *lcuth;
  TLine *lresi;

  int totalPRMC = h_decayPRMC->GetEntries();
  int totalNPMC = h_decayNPMC->GetEntries();

  for(int bins=0; bins<h_decayPRMC->GetNbinsX(); bins++){
    if(PRMC==0){
      h_deffPRMC->SetBinContent(bins,h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries());
      h_deffNPMC->SetBinContent(bins,1-(h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries()));

      if(h_decayPRMC->Integral(0,bins)<(0.899*totalPRMC)||h_decayPRMC->Integral(0,bins)>(totalPRMC*0.9099)) continue;
      cout<<bins<<"th bin"<<", ljpsi: "<<h_decayPRMC->GetBinCenter(bins)<<", PR eff: "<<h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries()*100<<"%"<<endl;
      lcutv = new TLine(h_decayPRMC->GetBinCenter(bins), 0, h_decayPRMC->GetBinCenter(bins), 1);
      lcuth = new TLine(xmin, h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries(), h_decayPRMC->GetBinCenter(bins), h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries());
      lresi = new TLine(xmin, h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries(), h_decayNPMC->GetBinCenter(bins), h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries());
    }
    else if(PRMC==1){
      h_deffNPMC->SetBinContent(bins,h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries());
      h_deffPRMC->SetBinContent(bins,1-(h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries()));

      if(h_decayPRMC->Integral(0,bins)<(0.97*totalPRMC)||h_decayPRMC->Integral(0,bins)>(totalPRMC*0.98)) continue;
      cout<<bins<<"th bin"<<", ljpsi: "<<h_decayPRMC->GetBinCenter(bins)<<", NP eff: "<<(1-h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries())*100<<"%"<<endl;
      lcutv = new TLine(h_decayPRMC->GetBinCenter(bins), 0, h_decayPRMC->GetBinCenter(bins), 1);
      lcuth = new TLine(xmin, h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries(), h_decayNPMC->GetBinCenter(bins), h_decayNPMC->Integral(0,bins)/h_decayNPMC->GetEntries());
      lresi = new TLine(xmin, h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries(), h_decayPRMC->GetBinCenter(bins), h_decayPRMC->Integral(0,bins)/h_decayPRMC->GetEntries());
    }
  }

  TCut ctauCut;
  if(PRMC==0){ctauCut = Form("ctau3D<=%.4f", lcutv->GetX1());}
  else if(PRMC==1){ctauCut = Form("ctau3D>=%.4f", lcutv->GetX1());}
  cout<<"LJpsi cut: "<<lcutv->GetX1()<<endl;
  treeData->Draw("mass>>h_mass",totalCut,"EP");
  treeData->Draw("mass>>h_massCut",totalCut&&ctauCut,"EP");
  cout<<"How many Jpsi w/o Cut??"<<h_mass->GetEntries()<<endl;
  cout<<"How many Jpsi w/ Cut??"<<h_massCut->GetEntries()<<endl;

  //gSystem->mkdir("Outputs/decayL",1);

  TCanvas* c_mass = new TCanvas("c_mass","",600,600);
  h_mass->Draw();
  h_mass->SetMinimum(0);
  h_mass->SetLineColor(kBlack);
  h_massCut->Draw("same");
  h_massCut->SetLineColor(kRed);
  c_mass->Update();
  if(PRMC==0){
  c_mass->SaveAs(Form("../figs/decayL/PRMC/mass_%s.pdf", kineLabel.Data()));}
  else if(PRMC==1){
  c_mass->SaveAs(Form("../figs/decayL/NPMC/mass_%s.pdf", kineLabel.Data()));}

  float pos_x = 0.23;
  float pos_x_mass = 0.55;
  float pos_y = 0.65;
  float pos_y_diff = 0.071;
  int text_color = 1;
  float text_size = 16;
  TString perc = "%";

  //TCanvas* c_decay_den = new TCanvas("c_decay_den","",600,600);
  //h_decay_den->Draw("P");
  ////lcutv->Draw("same");
  //c_decay_den->Update();
  //c_decay_den->SaveAs(Form("Outputs/makeV2Hist_%s/decayden_%s_%s.pdf",sample.Data(), consider.Data(), kineLabel.Data()));
  //if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  //else if(ptLow ==6.5) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  //else if(ptHigh==6.5) drawText(Form("%.f < p_{T}^{#mu#mu} < %.1f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  //else if(ptLow!=0) drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  //if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  //else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  //drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  //drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
  //drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);

  TCanvas* c_decayL = new TCanvas("c_decayL","",600,600);
  h_deffPRMC->Draw("l");
  h_deffNPMC->Draw("same");
  h_deffNPMC->SetLineColor(kRed+2);
  jumSun(xmin,1,xmax,1,1,1);
  lcutv->Draw("same");
  if(ptLow!=0) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff*0.5,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  //drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*1,text_color,text_size);
  drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*1.5,text_color,text_size);
  drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  if(PRMC==0){
    drawText(Form("L_{J/psi} cut: %.4f", lcutv->GetX1()),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
    drawText(Form("PR J/psi eff: %.3f", lcuth->GetY1()),pos_x_mass,pos_y-pos_y_diff*3.5,text_color,text_size);
    drawText(Form("NP J/psi res: %.3f", lresi->GetY1()),pos_x_mass,pos_y-pos_y_diff*4, kRed+2,text_size);
    c_decayL->Update();
    c_decayL->SaveAs(Form("../figs/decayL/PRMC/decay_%s.pdf",kineLabel.Data()));
    TFile *wf = new TFile(Form("../roots/decayL/PRMC/decay_hist_%s.root", kineLabel.Data()),"recreate");
    wf->cd();}
    else if(PRMC==1){
      drawText(Form("L_{J/psi} cut: %.4f", lcutv->GetX1()),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
      drawText(Form("NP J/psi eff: %.3f", 1-lcuth->GetY1()),pos_x_mass,pos_y-pos_y_diff*3.5,kRed+2,text_size);
      drawText(Form("PR J/psi res: %.3f", 1-lresi->GetY1()),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);
      c_decayL->Update();
      c_decayL->SaveAs(Form("../figs/decayL/NPMC/decay_%s.pdf",kineLabel.Data()));
      TFile *wf = new TFile(Form("../roots/decayL/NPMC/decay_hist_%s.root", kineLabel.Data()),"recreate");
      wf->cd();}

      h_deffPRMC->Write();
      h_deffNPMC->Write();
      h_mass->Write();
      h_massCut->Write();

      }
void GetHistSqrt(TH1D* h1, TH1D* h2){
  if(h1->GetNbinsX() != h2->GetNbinsX()){ cout << "Inconsistent # of bins b/w histograms !! " << endl;}
  double content;
  double err;
  for(int i=1; i<=h1->GetNbinsX(); i++){
    content=0;err=0;
    content = h1->GetBinContent(i);
    err = h1->GetBinError(i);
    err = 0.5*err*TMath::Power(content,-0.5);
    h2->SetBinContent(i,TMath::Sqrt(content));
    h2->SetBinError(i,err);
  }
} 
