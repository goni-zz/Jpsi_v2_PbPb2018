#include <iostream>
#include "rootFitHeaders.h"
#include "cutsAndBin.h"
#include "Style.h"
#include "commonUtility.h"
#include "tdrstyle.C"
#include "CMS_lumi_v2mass.C"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include <RooAbsData.h>
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "RooBinning.h"
#include "RooNumIntConfig.h"
#include "RooWorkspace.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "RooGaussModel.h"
#include "RooChi2Var.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooPlotable.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit;
void CtauErr()
{
  TFile *f1 = new TFile("MassFit/fitRes/fitresults_jpsi_DoubleCB_pt3.0-50.0_y1.6-2.4_muPt3.5_centrality20-180.root","READ");

  //	RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = (RooWorkspace*)f1->Get("workspace");
  //ws->import(*dataset);
  //	ws->data("dataset")->Print();
  cout << "J/psi yield : " << ws->var("nSig1s") -> getVal() << endl;
  cout << "Bkg yield : " << ws->var("nBkg") -> getVal() << endl;
  cout << "################################################" << endl;
  //	RooDataSet *reduceDS = (RooDataSet*)dataset->reduce(RooArgSet(*ws->var("ctau3DErr")) );
  //	ws->import(*reduceDS);
  //	ws->var("ctau3DErr")->setRange(0,0.2);

  RooRealVar *invMass=ws->var("mass");
  RooRealVar *sigYield = ws->var("nSig1s");
  RooRealVar *bkgYield = ws->var("nBkg");
  //RooDataSet* dataset = (RooDataSet*)ws->data("dataset")->Clone("TMP_DATA");
  RooDataSet* dataset = (RooDataSet*)ws->data("dataset");
  //RooAbsPdf* pdf = (RooAbsPdf*)(ws->pdf("model")->Clone("pdfMass"));
  RooAbsPdf* pdf = ws->pdf("model");
  RooRealVar *ctau3DErr = ws->var("ctau3DErr");

  RooArgList yieldList;
  yieldList.add(* ws->var("nSig1s"));
  yieldList.add(* ws->var("nBkg"));

  RooStats::SPlot sData = RooStats::SPlot("sData", "An SPlot", *dataset, pdf, yieldList);
  cout << "=========== SPlot done ===========" << endl;

  ws->import(*dataset,Rename("dataset_SPLOT"));

  cout << "J/psi yield : " << ws->var("nSig1s") -> getVal() << " , sWeights : " << sData.GetYieldFromSWeight("nSig1s") << endl;
  cout << "Bkg yield : " << ws->var("nBkg") -> getVal() << " , sWeights : " << sData.GetYieldFromSWeight("nBkg") << endl;

  ws->Print();
  dataset->Print();
  //	dataset_SPLOT->Print();
  //create weighted data sets
  double ctauErrMax = 1, ctauErrMin = -1;
  double binWidth = 0.0025;
  int nBins = 1000;

  TH1D* hTot = (TH1D*)ws->data("dataset")->createHistogram(("hTot"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrMin,ctauErrMax));
  hTot->Delete();
  cout << "HERE" << endl;

  RooDataSet* dataw_Bkg = new RooDataSet("TMP_BKG_DATA","TMP_BKG_DATA", RooArgSet(*ws->var("ctau3DErr"), *ws->var("nBkg_sw")));
  TH1D* hBkg = (TH1D*)ws->data("dataset")->createHistogram(("hBkg_sWeight"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrMin,ctauErrMax));

  RooDataSet* dataw_Jpsi = new RooDataSet("TMP_JPSI_DATA","TMP_JPSI_DATA", RooArgSet(*ws->var("ctau3DErr"), *ws->var("nSig1s_sw")));
  TH1D* hJpsi = (TH1D*)ws->data("dataset")->createHistogram(("hJpsi_sWeight"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrMin,ctauErrMax));

  TCanvas *cdata = new TCanvas("sPlot", "sPlot demo", 600, 600);
  //cdata->Divide(1, 3);
  //cdata->cd(1);
  //double minRange = (double)(floor(rangeErr[0]*100.)/100.);
  //double maxRange = (double)(ceil(rangeErr[1]*100.)/100.);
  RooPlot *frame = ctau3DErr->frame(Range(0, 0.22));
  dataset->plotOn(frame);
  //pdf->plotOn(frame, Name("Total pdf"));
  cdata->SetLogy();
  //pdf->plotOn(frame, Components(*zpdf), LineStyle(kDashed), LineColor(kRed), Name("Zpdf"));
  //pdf->plotOn(frame, Components(*qcdpdf), LineStyle(kDashed), LineColor(kGreen), Name("QCDpdf"));
  //TLegend leg(0.11, 0.5, 0.5, 0.8);
  //leg.AddEntry(frame->findObject("FullModel"), "Full model", "L");
  //leg.AddEntry(frame->findObject("ZModel"), "Z model", "L");
  //leg.AddEntry(frame->findObject("QCDModel"), "QCD model", "L");
  //leg.SetBorderSize(0);
  //leg.SetFillStyle(0);
  //frame->SetTitle("Fit of model to discriminating variable");
  //frame->SetMinimum(-1);
  //frame->SetMaximum(1);
  frame->Draw();
}

