#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../../commonUtility.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "../../Style.h"
#include "../../cutsAndBin.h"
#include "../../HiEvtPlaneList.h"
//#include "SimplePtFit.C"
using namespace std;

void vnlinearfit_Pt(
  float ptLow =  6.5, float ptHigh = 7.5,
  float yLow = 0, float yHigh = 2.4,
  int cLow = 40, int cHigh = 80,
  int weight_PR = 0, //PR : 0, NP : 1
  bool fEffW=true, bool fAccW=true, bool isPtW=true, bool isTnP=true
  ) {
  
  gStyle->SetOptStat(0);
  setTDRStyle();
  gStyle->SetOptFit(0000);
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;

  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0, cLow, cHigh);
  TString SimFit;
  if(cLow==20&&cHigh==120) SimFit = "10_60";
  else if(cLow==40&&cHigh==80) SimFit = "20_40";
  
  TFile *fFinal; 
  TFile *fv2fit1;
  TFile *fv2fit2;
  TFile *fv2fit3;
  
  fFinal = new TFile(Form("../../Macros/2021_09_14/roots/2DFit_%s/CtauRes_Sys/2DFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", SimFit.Data(), kineLabel.Data(), "PR", fEffW, fAccW, isPtW, isTnP));
  //if(ptLow==3&&ptHigh==6.5)fFinal = new TFile(Form("../../Macros/2021_09_14/roots/2DFit_%s/CtauRes_Sys/2DFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_v2.root", SimFit.Data(), kineLabel.Data(), "PR", fEffW, fAccW, isPtW, isTnP));
  //else fFinal = new TFile(Form("../../Macros/2021_09_14/roots/2DFit_%s/CtauRes_Sys/2DFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", SimFit.Data(), kineLabel.Data(), "PR", fEffW, fAccW, isPtW, isTnP));

  TH1D *Fraction1 = (TH1D*)fFinal->Get("Fraction1");
  TH1D *Fraction2 = (TH1D*)fFinal->Get("Fraction2");
  TH1D *Fraction3 = (TH1D*)fFinal->Get("Fraction3");

  double ctauLow =Fraction1->GetBinLowEdge((double)Fraction1->FindFirstBinAbove(1e-3))+Fraction1->GetBinWidth((double)Fraction1->FindFirstBinAbove(1e-3));
  double ctauHigh=Fraction2->GetBinLowEdge((double)Fraction2->FindFirstBinAbove(1e-3))+Fraction2->GetBinWidth((double)Fraction2->FindFirstBinAbove(1e-3));
  double bfrac1=Fraction1->GetBinContent((double)Fraction1->FindFirstBinAbove(1e-3));
  double bfrac2=Fraction2->GetBinContent((double)Fraction2->FindFirstBinAbove(1e-3));
  double bfrac3=Fraction3->GetBinContent((double)Fraction3->FindFirstBinAbove(1e-3));
  
  fv2fit1=new TFile(Form("roots/SignalPDFvariation/%s/SimFitResult_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.root", SimFit.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLow,"ctauL"));
  fv2fit2=new TFile(Form("roots/SignalPDFvariation/%s/SimFitResult_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%.5f_%s.root", SimFit.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLow,ctauHigh,"ctauC"));
  fv2fit3=new TFile(Form("roots/SignalPDFvariation/%s/SimFitResult_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.root", SimFit.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauHigh,"ctauR"));
  //fv2fit3=new TFile(Form("../simultaneous_fit/fitting/roots/v2mass_fit/20_40/SimFitResult_pt3.0-4.5_y1.6-2.4_muPt0.0_centrality40-80_Eff1_Acc1_PtW1_TnP1_ctau_0.13000_ctauR.root"));

  //fv2fit1=new TFile(Form("../simultaneous_fit/fitting/roots/v3mass_fit/%s/SimFitResult_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.root", SimFit.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLow,"ctauL"));
  //fv2fit2=new TFile(Form("../simultaneous_fit/fitting/roots/v3mass_fit/%s/SimFitResult_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%.5f_%s.root", SimFit.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLow,ctauHigh,"ctauC"));
  //fv2fit3=new TFile(Form("../simultaneous_fit/fitting/roots/v3mass_fit/%s/SimFitResult_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.root", SimFit.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauHigh,"ctauR"));
  double v2BinP1,v2BinP2,v2BinP3;
  double v2BinE1,v2BinE2,v2BinE3;
  double ptBin=(ptLow+ptHigh)/2;
  TGraphErrors *v2vspt1 = (TGraphErrors*)fv2fit1->Get("v2vspt");
  TGraphErrors *v2vspt2 = (TGraphErrors*)fv2fit2->Get("v2vspt");
  TGraphErrors *v2vspt3 = (TGraphErrors*)fv2fit3->Get("v2vspt");
  v2vspt1->GetPoint(0,ptBin,v2BinP1);
  v2vspt2->GetPoint(0,ptBin,v2BinP2);
  v2vspt3->GetPoint(0,ptBin,v2BinP3);
  v2BinE1=v2vspt1->GetErrorY(0);
  v2BinE2=v2vspt2->GetErrorY(0);
  v2BinE3=v2vspt3->GetErrorY(0);

  double bfracBin[] = {bfrac1, bfrac3, bfrac2};
  double bfracErr[] = {0.0, 0.0, 0.0};
  double v2Bin[] = {v2BinP1,v2BinP2,v2BinP3};
  double v2Err[] = {v2BinE1,v2BinE2,v2BinE3};

  TGraphErrors *gFrac_v2 = new TGraphErrors(3, bfracBin, v2Bin, bfracErr, v2Err);
  TF1 *fit1 = new TF1("fit1", "[0]*x +[1]", -10.,10.);
  fit1->SetLineStyle(1);
  fit1->SetLineColor(kRed+2);
  fit1->SetLineWidth(3);
  gFrac_v2->SetMarkerStyle(20);
  gFrac_v2->SetMarkerColor(kBlue+2);
  gFrac_v2->SetLineColor(kBlue+2);
  gFrac_v2->SetMarkerSize(1.4);
  gStyle->SetPadTickY(1);
  gFrac_v2->Fit("fit1", "R", "R", 0., 1.);
  ///gFrac_v2->GetYaxis()->SetRangeUser(0, 0.1);
  //handsomeTH1(gFrac_v2,1);
  TCanvas*c1= new TCanvas("c1","c1",1000,800);
  TPad* pad1 = new TPad("pad1","pad1",0,0,1.0,1.0);
  c1->cd();
  pad1->SetTicks(1,1);
  //pad1->SetBottomMargin(0);
  //pad1->SetLeftMargin(0.19);
  //pad1->SetTopMargin(0.08);
  pad1->Draw();
  pad1->cd();
  double pad1W = pad1->GetWw()*pad1->GetAbsWNDC();
  double pad1H = pad1->GetWh()*pad1->GetAbsHNDC();
  double tickScaleX = (pad1->GetUxmax() - pad1->GetUxmin())/(pad1->GetX2()-pad1->GetX1())*pad1H;
  TH1D *hint = new TH1D("hint", "", 1000, 0., 1.);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint, 0.68);
  //gFrac_v2->Draw("p");
  hint->SetStats(kTRUE);
  hint->SetFillColor(kBlue-10);
  hint->SetFillStyle(1001);
  hint->Draw("e3 same");
  hint->GetYaxis()->SetRangeUser(0,0.2);
  hint->GetYaxis()->SetTitle("v_{2}^{Sig}");
  hint->GetYaxis()->SetTitleSize(0.05);
  hint->GetYaxis()->SetLabelSize(0.04);
  hint->GetYaxis()->CenterTitle();
  hint->GetYaxis()->CenterTitle();
  //hint->GetYaxis()->SetTitleFont(42);
  hint->GetXaxis()->SetTitle("Nonprompt J/#psi fraction");
  hint->GetXaxis()->SetLabelSize(0.04);
  hint->GetXaxis()->SetTitleSize(0.05);
  hint->GetXaxis()->CenterTitle();
  //hint->GetXaxis()->SetLabelFont(42);
  gFrac_v2->Draw("p same");
  pad1->Update();
  CMS_lumi_v2mass(pad1,iPeriod,iPos);

  cout<<" prompt J/psi v2 : "<<fit1->Eval(0.0)<<" +- "<<hint->GetBinError(1)<<endl;
  cout<<" Nonprompt J/psi v2 : "<<fit1->Eval(1.0)<<" +- "<<hint->GetBinError(1000)<<endl;
  double v2prp = (double)fit1->Eval(0.0);
  double v2eprp = (double)hint->GetBinError(1);
  double v2nprp = (double)fit1->Eval(1.0);
  double v2enprp = (double)hint->GetBinError(1000);
  //drawText(Form("3 < p_{T}^{#mu#mu} < 4.5 GeV/c"),0.73,0.76,1,24);
  //drawText(Form("4.5 < p_{T}^{#mu#mu} < 6.5 GeV/c"),0.73,0.76,1,24);
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow,ptHigh),0.73,0.76,1,24);
  //drawText(Form("7.5 < p_{T}^{#mu#mu} < 9 GeV/c"),0.73,0.76,1,24);
  //drawText(Form("9 < p_{T}^{#mu#mu} < 12 GeV/c"),0.73,0.76,1,24);
  //drawText(Form("12 < p_{T}^{#mu#mu} < 50 GeV/c"),0.73,0.76,1,24);
  drawText(Form("|y^{#mu#mu}| < 2.4"),0.73,0.70,1,24);
  drawText(Form("Centrality %i-%i%%",cLow/2,cHigh/2),0.73,0.64,1,24);
  drawText(Form("Prompt J/#psi v_{2} = %0.3f #pm %0.3f",v2prp,v2eprp),0.23,0.26,1,26);
  drawText(Form("Non-Prompt J/#psi v_{2} = %0.3f #pm %0.3f",v2nprp,v2enprp),0.23,0.20,1,26);

  c1->SaveAs(Form("SignalPDFvariation_fraction_vs_v2_%s.pdf",kineLabel.Data()));
  TH1D *prv2=new TH1D("prv2","",1,0,1);
  prv2->SetBinContent(0,v2prp);
  prv2->SetBinError(0,v2eprp);
  TH1D *npv2=new TH1D("npv2","",1,0,1);
  npv2->SetBinContent(0,v2nprp);
  npv2->SetBinError(0,v2enprp);
  TFile *fv2 = new TFile(Form("SignalPDFvariation_v2_%s.root",kineLabel.Data()),"RECREATE");
  fv2->cd();
  prv2->Write();
  npv2->Write();
}
