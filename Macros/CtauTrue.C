#include <iostream>
#include "../rootFitHeaders.h"
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBin.h"
#include "../CMS_lumi_v2mass.C"
#include "../tdrstyle.C"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit;

void CtauTrue(
    float ptLow=8, float ptHigh=12,
    float yLow=0, float yHigh=2.4,
    int cLow=20, int cHigh=120,
    float muPtCut=0.0,
    bool whichModel=0,  // Nominal=0, Alternative=1
    int ICset=1,
    int PR=2, //0=PR, 1=NP, 2=Inc.
    float ctauCut=0.08
    )
{
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir("../roots/2DFit_1203/");
  gSystem->mkdir("../figs/2DFit_1203/");

  TString bCont;
  if(PR==0) bCont="Prompt";
  else if(PR==1) bCont="NonPrompt";
  else if(PR==2) bCont="Inclusive";

  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  TFile* f1; TFile* f2; 
  TString kineCut;
  TString kineCutMC;
  TString SigCut;
  TString BkgCut;
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

  f1 = new TFile("../skimmedFiles/OniaRooDataSet_NonPrompt_GenInReco_3_6p5_201005.root");
  f2 = new TFile("../skimmedFiles/OniaRooDataSet_NonPrompt_GenInReco_6p5_50_201005.root");
  if(PR==2) {
    kineCutMC = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.5" ,ptLow, ptHigh, yLow, yHigh);
  }
  else if(PR==0) {
    kineCutMC = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f &&  mass>2.6 && mass<3.5 && ctau3Dtrue<%.2f",ptLow, ptHigh, yLow, yHigh, ctauCut);
  }
  // SigCut  = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && cBin>%d && cBin<%d",ptLow, ptHigh, yLow, yHigh, cLow, cHigh);
  //BkgCut  = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && ((mass>2.6 && mass <= 2.8) || (mass>=3.2&&mass<3.5)) && cBin>%d && cBin<%d",ptLow, ptHigh, yLow, yHigh, cLow, cHigh);
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

  kineCutMC = accCut+kineCutMC;
  //SigCut = SigCut;
  //BkgCut = BkgCut;

  //RooDataSet *datasetMC = (RooDataSet*)f2->Get("dataset");
  RooDataSet *datasetMC;
  if(ptHigh <= 6.5) {datasetMC = (RooDataSet*)f1->Get("dataset");}
  else {datasetMC = (RooDataSet*)f2->Get("dataset");}

  RooWorkspace *ws = new RooWorkspace("workspace");
  RooWorkspace *wsmc = new RooWorkspace("workspaceMC");
  wsmc->import(*datasetMC);
  //wsmc->import(*datasetMC2);
  RooCategory tp("tp","tp");
  tp.defineType("C");
  //tp.defineType("D");
  RooDataSet *reducedDS_MC = (RooDataSet*)datasetMC->reduce(RooArgSet(*(wsmc->var("ctau3Dtrue")), *(wsmc->var("mass")), *(wsmc->var("pt")), *(wsmc->var("y"))), kineCutMC.Data() );
  //RooDataSet* reducedDS_MC = new RooDataSet("reducedDS_MC","reducedDS_MC",RooArgSet(*(wsmc->var("ctau3Dtrue")), *(wsmc->var("mass")), *(wsmc->var("pt")), *(wsmc->var("y"))),Index(tp),Import("C",*datasetMC1),Import("D",*datasetMC2));
  reducedDS_MC->SetName("reducedDS_MC");
  reducedDS_MC->Print();
  ws->import(*reducedDS_MC);
  //ws->var("ctau3Dtrue")->setRange(0, 5);
  ws->var("ctau3Dtrue")->setRange("ctauTrueRange", 0, 6.);
  ws->var("ctau3Dtrue")->Print();

  //***********************************************************************
  //**************************** CTAU TRUE FIT ****************************
  //***********************************************************************
  cout << endl << "************ Start MC Ctau True Fit ***************" << endl << endl;
  //MC NP ctau true
  double entries_True = ws->data("reducedDS_MC")->numEntries();
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, entries_True*0.9, entries_True));
  ws->factory("lambdaDSS[0.1, 1e-6, 10.]");
  // create the PDF
  ws->factory(Form("TruthModel::%s(%s)", "TruthModel_ctauTrue", "ctau3Dtrue"));
  //ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "ctauTruePdf", "ctau3Dtrue",  
  //      "lambdaDSS", 
  //      "TruthModel_ctauTrue"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "ctauTruePdf",
        "ctau3Dtrue",
        "lambdaDSS",
        "TruthModel_ctauTrue"));

  RooAbsPdf *ctauTrueModel = ctauTrueModel = new RooAddPdf("TrueModel_Tot", "TrueModel_Tot", *ws->pdf("ctauTruePdf"), *ws->var("N_Jpsi_MC"));
  ws->import(*ctauTrueModel);
  //TH1 *hCTrue = (TH1*)ctauTrueModel->createHistogram("ctau3Dtrue",50,50);

  TCanvas* c_D =  new TCanvas("canvas_D","My plots",4,565,550,520);
  c_D->cd();
  TPad *pad_D_1 = new TPad("pad_D_1", "pad_D_1", 0, 0.16, 0.98, 1.0);
  pad_D_1->SetTicks(1,1);
  pad_D_1->Draw(); pad_D_1->cd();
  RooPlot* myPlot_D = wsmc->var("ctau3Dtrue")->frame(Bins(nCtauTrueBins), Range(-1, 6)); // bins
  ws->data("reducedDS_MC")->plotOn(myPlot_D,Name("mcHist_D"));
  myPlot_D->SetTitle("");


  pad_D_1->cd();
  gPad->SetLogy();
  RooFitResult* fitCtauTrue = ws->pdf("TrueModel_Tot")->fitTo(*reducedDS_MC, SumW2Error(false), Range("-1e-6, 6"), Extended(kTRUE), NumCPU(nCPU), PrintLevel(3), Save());
  RooPlot* myPlot2_D = (RooPlot*)myPlot_D->Clone();
  ws->data("reducedDS_MC")->plotOn(myPlot2_D,Name("MCHist_Tot"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), Binning(nCtauTrueBins));
  ws->pdf("TrueModel_Tot")->plotOn(myPlot2_D,Name("MCpdf_Tot"), Normalization(ws->data("reducedDS_MC")->sumEntries(), RooAbsReal::NumEvent),
      Precision(1e-4), LineColor(kRed+2), NormRange("ctauTrueRange"), Range("ctauTrueRange"));
  myPlot2_D->GetYaxis()->SetRangeUser(10e-2, ws->data("reducedDS_MC")->sumEntries()*10);
  myPlot2_D->GetXaxis()->SetRangeUser(-1, 7);
  myPlot2_D->GetXaxis()->CenterTitle();
  myPlot2_D->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} MC True (mm)");
  myPlot2_D->SetFillStyle(4000);
  myPlot2_D->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_D->GetXaxis()->SetLabelSize(0);
  myPlot2_D->GetXaxis()->SetTitleSize(0);
  myPlot2_D->Draw();

  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x+0.5,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x+0.5,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x+0.5,text_y-y_diff,text_color,text_size);
  drawText(Form("#lambdaDSS = %.4f #pm %.4f", ws->var("lambdaDSS")->getVal(), ws->var("lambdaDSS")->getError() ),text_x+0.5,text_y-y_diff*2,text_color,text_size);

  TPad *pad_D_2 = new TPad("pad_D_2", "pad_D_2", 0, 0.006, 0.98, 0.227);
  RooPlot* frameTMP_D = (RooPlot*)myPlot2_D->Clone("TMP");
  RooHist* hpull_D;
  pullDist(ws, pad_D_2, c_D, frameTMP_D, hpull_D, "MCHist_Tot", "MCpdf_Tot", "ctau3Dtrue", nCtauTrueBins, -1, 6, "#font[12]{l}_{J/#psi} MC True (mm)");
  printChi2_test(ws, pad_D_2, frameTMP_D, hpull_D, fitCtauTrue, "ctau3Dtrue", "MCHist_Tot", "MCpdf_Tot", nCtauTrueBins);

  TLine *lD = new TLine(-1, 0, 6.5, 0);
  lD->SetLineStyle(1);
  lD->Draw("same");
  //printChi2(ws, pad_D_2, frameTMP_D, fitCtauTrue, "ctau3Dtrue", "MCHist_Tot", "MCpdf_Tot", nCtauTrueBins, false);
  pad_D_2->Update();

  c_D->Update();
  c_D->SaveAs(Form("../figs/2DFit_1203/ctauTrue_%s_%s.pdf",bCont.Data(),kineLabel.Data()));

  //TH1 *h1 = (TH1*)TrueModel_Tot->createHistogram("ctau3Dtrue",50,50);
  TFile *outFile = new TFile(Form("../roots/2DFit_1203/CtauTrueResult_%s_%s.root",bCont.Data(),kineLabel.Data()),"RECREATE");
  RooArgSet* fitargs = new RooArgSet();
  fitargs->add(fitCtauTrue->floatParsFinal());
  ctauTrueModel->Write();
  outFile->Close();
  //TFile *hFile = new TFile(Form("ctauTreuHist_%s.root",kineLabel.Data()),"recreate");
  //hCTrue->Write();
  //hFile->Close();
  cout << endl << "************ Finished MC Ctau True Fit ***************" << endl << endl;

}
