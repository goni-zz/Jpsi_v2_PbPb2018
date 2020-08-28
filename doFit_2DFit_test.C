//This code fits the Jpsi data with either the nominal fit or an alternative fit. The difference between the two fits is the signal shape. The nominal fit fits the signals with double CB functions, while the alternative fit fits them with just a gaussian.

#include <iostream>
#include "rootFitHeaders.h"
#include "commonUtility.h"
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
#include "cutsAndBin.h"
#include "CMS_lumi_v2mass.C"
#include "tdrstyle.C"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit;
void doFit_2DFit_test( 
       float ptLow=6.5, float ptHigh=8.5, 
       float yLow=1.6, float yHigh=2.4,
       int cLow=20, int cHigh=120,
       float muPtCut=0.0,
       bool whichModel=0,   // Nominal = 0. Alternative = 1.
       int ICset = 1
			) 
{
  gStyle->SetEndErrorSize(0);

  float massLow = 2.6, massHigh = 3.5;
  int   nMassBin  = 36; //(massHigh-massLow)*30;

  float ctauLow = -3, ctauHigh = 5;
  float ctauResLow = -20, ctauResHigh = 20;
  int   nCtauBins  = (ctauHigh-ctauLow)*10;
  int   nCtauErrBins = 72;
  int   nCtauResBins = 72;
  //int nCtauErrBins = 100;
  
  double ctauErrLow = 1e-6, ctauErrHigh = 0.18;
  double binWidth = 0.0025;

  TFile* f1; TFile* f2;
  TString kineCut;
  TString SigCut;
  TString BkgCut;

  //Select Data Set
    //f1 = new TFile("skimmedFiles/OniaRooDataSet_isMC0_JPsi.root");
    f1 = new TFile("skimmedFiles/OniaRooDataSet_isMC0_JPsi1SW_2020819.root");
    f2 = new TFile("skimmedFiles/OniaRooDataSet_isMC1_JPsi1SW_2020821.root");
    kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.5",ptLow, ptHigh, yLow, yHigh);
    SigCut  = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.8 && mass<3.2",ptLow, ptHigh, yLow, yHigh);
    BkgCut  = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && ((mass>2.6 && mass <= 2.8) || (mass>=3.2&&mass<3.5))",ptLow, ptHigh, yLow, yHigh);
    //SigCut = Form("pt>%.2f && pt<%.2f && abs(y)<%.2f && abs(eta1)<%.2f && abs(eta2)<%.2f && mass>2.9&&mass<3.2", ptLow, ptHigh, yHigh, yHigh, yHigh);
    //BkgCut = Form("pt>%.2f && pt<%.2f && abs(y)<%.2f && abs(eta1)<%.2f && abs(eta2)<%.2f && (mass<2.9||mass>3.2)",ptLow, ptHigh, yHigh, yHigh, yHigh);
    TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

  //kineCut = accCut + kineCut;
  //SigCut = accCut + SigCut;
  //BkgCut = accCut + BkgCut;
  
  kineCut = kineCut;
  SigCut = SigCut;
  BkgCut = BkgCut;
  //import and merge datasets
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooDataSet *datasetMC = (RooDataSet*)f2->Get("dataset");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  cout << "####################################" << endl;
  RooDataSet *reducedDS_A = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), SigCut.Data() );
  RooDataSet *reducedDS_B = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), BkgCut.Data() );
  reducedDS_A->SetName("reducedDS_A");
  reducedDS_B->SetName("reducedDS_B");
  //ws->import(*reducedDS_A);
  //ws->import(*reducedDS_B);
  //ws->var("mass")->setRange(massLow, massHigh);
  //ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  //ws->var("ctau3DErr")->setRange(ctauLow, ctauHigh);
  //ws->var("ctau3DRes")->setRange(ctauResLow, ctauResHigh);
  //reducedDS_A->Print();
  //reducedDS_B->Print();
  RooCategory tp("tp","tp");
  tp.defineType("A");
  tp.defineType("B");

  // Create a dataset that imports contents of all the above datasets mapped by index category tp
  //RooDataSet* dsAB = new RooDataSet("dsAB","dsAB",RooArgSet(*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))),Index(tp),Import("A",*reducedDS_A),Import("B",*reducedDS_B));
  RooDataSet* dsAB = new RooDataSet("dsAB","dsAB",RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))),Index(tp),Import("A",*reducedDS_A),Import("B",*reducedDS_B));
  cout << "******** New Combined Dataset ***********" << endl;
  dsAB->Print();
  ws->import(*dsAB);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  ws->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);
  ws->var("ctau3DErr")->setRange(ctauErrLow, ctauErrHigh);
  ws->var("ctau3DErr")->setRange("ctauErrRange",ctauErrLow, ctauErrHigh);
  ws->var("ctau3DRes")->setRange(ctauResLow, ctauResHigh);
  ws->var("ctau3DRes")->setRange("ctauResRange", ctauResLow, ctauResHigh);
  ws->var("mass")->Print();
  ws->var("ctau3D")->Print();
  ws->var("ctau3DErr")->Print();

  RooWorkspace *wsmc = new RooWorkspace("workspaceMC");
  wsmc->import(*datasetMC);
  RooDataSet *reducedDS_MC = (RooDataSet*)datasetMC->reduce(RooArgSet(*(wsmc->var("ctau3Dtrue")), *(wsmc->var("mass")), *(wsmc->var("pt")), *(wsmc->var("y"))), kineCut.Data() );
  reducedDS_MC->SetName("reducedDS_MC");
  wsmc->import(*reducedDS_MC);
  reducedDS_MC->Print();
  wsmc->var("ctau3Dtrue")->setRange(0, 5);
  wsmc->var("ctau3Dtrue")->setRange("ctauTrueRange", 0, 5);
  wsmc->var("ctau3Dtrue")->Print();
  
  TCanvas* c_A =  new TCanvas("canvas_A","My plots",4,45,550,520);
  c_A->cd();
  RooPlot* myPlot_A = ws->var("mass")->frame(nMassBin); // bins
  ws->data("dsAB")->plotOn(myPlot_A,Name("dataHist_A"));

  TCanvas* c_B =  new TCanvas("canvas_B","My plots",504,45,550,520);
  c_B->cd();
  RooPlot* myPlot_B = ws->var("ctau3DErr")->frame(nCtauErrBins); // bins
 // ws->data("dsAB")->plotOn(myPlot_B,Name("dataHist_B"));

  TCanvas* c_C =  new TCanvas("canvas_C","My plots",504,45,550,520);
  c_C->cd();
  RooPlot* myPlot_C = ws->var("ctau3DRes")->frame(nCtauResBins); // bins
  //ws->data("dsAB")->plotOn(myPlot_C,Name("dataHist_C"));

  TCanvas* c_D =  new TCanvas("canvas_D","My plots",504,45,550,520);
  c_D->cd();
  RooPlot* myPlot_D = wsmc->var("ctau3Dtrue")->frame(nCtauBins); // bins
  wsmc->data("reducedDS_MC")->plotOn(myPlot_D,Name("mcHist_D"));

  TCanvas* c_E =  new TCanvas("canvas_E","My plots",504,45,550,520);
  c_E->cd();
  RooPlot* myPlot_E = ws->var("ctau3D")->frame(nCtauBins); // bins
  //ws->data("dsAB")->plotOn(myPlot_C,Name("dataHist_C"));

  TCanvas* c_F =  new TCanvas("canvas_F","My plots",4,45,550,520);
  c_F->cd();
  RooPlot* myPlot_F = ws->var("ctau3D")->frame(nCtauBins); // bins
  //ws->data("dsAB")->plotOn(myPlot_F,Name("dataHist_ctauTot"));

  //***********************************************************************
  //****************************** MASS FIT *******************************
  //***********************************************************************
  //         The order is {sigma_1,  x, alpha_1, n_1,   f, err_mu, err_sigma, m_lambda}
  double paramsupper[8] = {0.2,    3.0,   3.321, 5.0, 1.0,   25.0,      25.0,     25.0};
  double paramslower[8] = {0.02,   0.0,     1.0, 1.0, 0.0,    0.0,       0.0,      0.0};
  //SIGNAL: initial params
  double sigma_1_init = 0.1;
  double x_init = 0.8;
  double alpha_1_init = 1.5;
  double n_1_init = 2.6;
  double f_init = 0.5;
  if (ICset>1 && ICset<4) {
    sigma_1_init = 0.3;
    x_init = 0.3;
    alpha_1_init = 2.6;
    n_1_init = 3.0;
    f_init = 0.1;
  }

//SIGNAL
  RooRealVar    mean("m_{J/#Psi}","mean of the signal gaussian mass PDF",pdgMass.JPsi, pdgMass.JPsi -0.1, pdgMass.JPsi + 0.1 ) ;
  RooRealVar    sigma_1_A("sigma_1_A","width/sigma of the signal gaussian mass PDF",sigma_1_init, paramslower[0], paramsupper[0]);
  RooRealVar   *x_A = new RooRealVar("x_A","sigma ratio ", x_init, paramslower[1], paramsupper[1]);
  RooFormulaVar sigma_2_A("sigma_2_A","@0*@1",RooArgList(sigma_1_A, *x_A) );
  RooRealVar    alpha_1_A("alpha_1_A","tail shift", alpha_1_init, paramslower[2], paramsupper[2]);
  RooFormulaVar alpha_2_A("alpha_2_A","1.0*@0",RooArgList(alpha_1_A) );
  RooRealVar    n_1_A("n_1_A","power order", n_1_init , paramslower[3], paramsupper[3]);
  RooFormulaVar n_2_A("n_2_A","1.0*@0",RooArgList(n_1_A) );
  RooRealVar   *f = new RooRealVar("f","cb fraction", f_init, paramslower[4], paramsupper[4]);
  //Set up crystal ball shapes
  RooCBShape* cb_1_A = new RooCBShape("cball_1_A", "cystal Ball", *(ws->var("mass")), mean, sigma_1_A, alpha_1_A, n_1_A);
  RooAddPdf*  cb_A;
  //DOUBLE CRYSTAL BALL
  RooCBShape* cb_2_A = new RooCBShape("cball_2_A", "cystal Ball", *(ws->var("mass")), mean, sigma_2_A, alpha_2_A, n_2_A);
  cb_A = new RooAddPdf("cb_A","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );
  RooRealVar *nSig= new RooRealVar("nSig","inclusive Jpsi signals",0,100000);
  //RooRealVar *nSigPR= new RooRealVar("nSigPR","prompt Jpsi signals",0,100000);
  //RooRealVar *nSigNP= new RooRealVar("nSigNP","non-prompt Jpsi signals",0,100000);
  cb_A = new RooAddPdf("cb_A","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );
//BACKGROUND
  double err_sigma_init = 5;
  double err_mu_init = 8;
  double m_lambda_init = 5;
  if (ICset>2) {
    err_mu_init = 5;
    m_lambda_init = 5;
  }

//A: Mass Bkg function
  RooRealVar err_mu_A("#mu_A","err_mu", err_mu_init,  paramslower[5], paramsupper[5]) ;
  RooRealVar err_sigma_A("#sigma_A","err_sigma", err_sigma_init, paramslower[6], paramsupper[6]);
  RooRealVar m_lambda_A("#lambda_A","m_lambda",  m_lambda_init, paramslower[7], paramsupper[7]);
  //THIS IS THE BACKGROUND FUNCTION
  RooGenericPdf *bkg = new RooGenericPdf("bkg","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_A));
  RooRealVar *nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",10000,0,50000); 

  //RooRealVar *sl1 = new RooRealVar("sl1","sl1",0.1,-3.,3);
  //RooRealVar *cnst1 = new RooRealVar("cnst1","cnst1",0.1,-2.,2);
  //RooGenericPdf *bkg_1order = new RooGenericPdf("bkg","Background","@0*@1+@2",RooArgList( *(ws->var("mass")), sl1, cnst1) );
  //RooChebychev *bkg_1order;
  //bkg_1order = new RooChebychev("bkg","Background",*(ws->var("mass")),RooArgList(*sl1,*cnst1));
  //Build the model
//Model A: Mass
  RooAddPdf* model_A = new RooAddPdf();
  model_A = new RooAddPdf("model_A","Jpsi + Bkg",RooArgList(*cb_A, *bkg),RooArgList(*nSig,*nBkg));
  //model_A = new RooAddPdf("model_A","Jpsi + Bkg",RooArgList(*cb_A, *bkg_1order),RooArgList(*nSig,*nBkg));
  //model_A = new RooAddPdf("model_A","PR Jpsi + NP Jpsi + Bkg",RooArgList(*cb_1_A, *cb_2_A, *bkg),RooArgList(*nSigPR,*nSigNP,*nBkg));
  ws->import(*model_A);
  cout << endl << "********* Starting Mass Dist. Fit **************" << endl << endl;
  RooFitResult* fitResSim = ws->pdf("model_A")->fitTo(*dsAB,Save(), Hesse(kTRUE), Range(massLow,massHigh), Timer(kTRUE), Extended(kTRUE), NumCPU(12), PrintLevel(-1));
  cout << endl << "********* Finished Mass Dist. Fit **************" << endl << endl;

  c_A->cd();
  gPad->SetLogy();
  RooPlot* myPlot2_A = (RooPlot*)myPlot_A->Clone();
  dsAB->plotOn(myPlot2_A,Name("dataOS_FIT_A"),MarkerSize(.8));
  ws->pdf("model_A")->plotOn(myPlot2_A,Name("modelHist_A"), LineColor(kBlack));
  //ws->pdf("model_A")->plotOn(myPlot2_A,Name("Sig_A"),Components(RooArgSet(*cb_A)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_A")->plotOn(myPlot2_A,Name("bkgPDF_A"),Components(RooArgSet(*bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  myPlot2_A->SetFillStyle(4000);
  myPlot2_A->GetYaxis()->SetTitleOffset(1.43);
  //myPlot2_A->GetYaxis()->CenterTitle();
  //myPlot2_A->GetYaxis()->SetTitleSize(0.058);
  //myPlot2_A->GetYaxis()->SetLabelSize(0.054);
  //myPlot2_A->GetYaxis()->SetRangeUser(2*10, 2*10e3);
  myPlot2_A->GetYaxis()->SetRangeUser(ws->var("nSig")->getVal()/100, ws->var("nSig")->getVal());
  //myPlot2_A->SetMinimum(2*10);
  //myPlot2_A->GetXaxis()->SetLabelSize(0);
  myPlot2_A->GetXaxis()->CenterTitle();
  myPlot2_A->GetXaxis()->SetRangeUser(massLow,massHigh);
  myPlot2_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  //myPlot2_A->GetXaxis()->SetTitleSize(0);
  myPlot2_A->Draw();

  float text_x = 0.15;
  float text_y = 0.816;
  float y_diff = 0.05;
  float text_size = 15;
  int text_color = 1;
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);
  drawText(Form("n_{J/#psi} = %.f #pm %.f",ws->var("nSig")->getVal(),ws->var("nSig")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
  drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("nBkg")->getVal(),ws->var("nBkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
  //***********************************************************************
  //**************************** CTAU ERR FIT *****************************
  //***********************************************************************
  cout << endl << "*************** Start SPLOT *****************" << endl << endl;
  //SPlot Ctau Error
  RooRealVar *sigYield = ws->var("nSig");
  RooRealVar *bkgYield = ws->var("nBkg");
  RooArgList yieldList;
  yieldList.add(*ws->var("nSig"));
  yieldList.add(*ws->var("nBkg"));
  cout<<"Sig Yield: "<<sigYield->getVal()<<" +/- "<<sigYield->getError()<<endl;
  cout<<"Bkg Yield: "<<bkgYield->getVal()<<" +/- "<<bkgYield->getError()<<endl;

  c_B->cd();
  c_B->SetLogy();
  ////RooPlot *frame = ws->var("ctau3DErr")->frame(Range(0, ctauErrHigh)); //Bins(100);
  //original
  /*RooPlot* myPlot2_B = (RooPlot*)myPlot_B->Clone();
  TH1D* hTot = (TH1D*)ws->data("dsAB")->createHistogram(("hTot"), *ws->var("ctau3DErr"),Binning(nCtauErrBins,ctauErrLow,ctauErrHigh));
  TH1D* hSig = (TH1D*)ws->data("reducedDS_A")->createHistogram(("hSig"), *ws->var("ctau3DErr"),Binning(nCtauErrBins,ctauErrLow,ctauErrHigh));
  TH1D* hBkg = (TH1D*)ws->data("reducedDS_B")->createHistogram(("hBkg"), *ws->var("ctau3DErr"),Binning(nCtauErrBins,ctauErrLow,ctauErrHigh));
  RooDataHist* totHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hTot);
  RooDataHist* sigHist = new RooDataHist("reducedDS_A", "", RooArgSet(*ws->var("ctau3DErr")), hSig);
  RooDataHist* bkgHist = new RooDataHist("reducedDS_B", "", RooArgSet(*ws->var("ctau3DErr")), hBkg);
  RooHistPdf* TotPdf = new RooHistPdf("TotPdf","hist pdf", *ws->var("ctau3DErr"), *totHist);
  RooHistPdf* sigPdf = new RooHistPdf("sigPdf","hist pdf", *ws->var("ctau3DErr"), *sigHist);
  RooHistPdf* bkgPdf = new RooHistPdf("bkgPdf","hist pdf", *ws->var("ctau3DErr"), *bkgHist);
  RooAddPdf* model_ctErr = new RooAddPdf();
  model_ctErr = new RooAddPdf("model_ctErr","Sig + Bkg ctau models", RooArgList(*sigPdf, *bkgPdf), RooArgList(*sigYield, *bkgYield));
  ws->import(*model_ctErr);
  RooDataSet* data = (RooDataSet*)ws->data("dsAB")->Clone("TMP_DATA");
  RooArgSet* cloneSet = (RooArgSet*)RooArgSet(*ws->pdf("model_A"),"model_A").snapshot(kTRUE);
  auto pdf = (RooAbsPdf*)cloneSet->find("model_A");
  pdf->setOperMode(RooAbsArg::ADirty, kTRUE);
  RooStats::SPlot sData = RooStats::SPlot("sData", "An SPlot", *data, pdf, RooArgList(*sigYield, *bkgYield));
  ws->import(*dsAB, Rename("dataset_SPLOT"));

  cout<<"#### nSig: "<<ws->var("nSig")->getVal()<<", sWeights :"<<sData.GetYieldFromSWeight("nSig")<<endl;
  cout<<"#### nBkg: "<<ws->var("nBkg")->getVal()<<", sWeights :"<<sData.GetYieldFromSWeight("nBkg")<<endl;

  ws->data("dsAB")->plotOn(myPlot2_B,Name("dataHist_Tot"), MarkerSize(.7), Binning(nCtauErrBins));
  ws->pdf("model_ctErr")->plotOn(myPlot2_B,Name("totPdf"), LineColor(kGreen+1), Range("ctauErrRange"), LineWidth(2));
  ws->data("reducedDS_A")->plotOn(myPlot2_B,Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2), Binning(nCtauErrBins));
  ws->pdf("sigPdf")->plotOn(myPlot2_B,Name("sigPdf"),LineColor(kRed+2), LineWidth(2), Range("ctauErrRange"));
  ws->data("reducedDS_B")->plotOn(myPlot2_B,Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue+2), MarkerColor(kBlue+2), Binning(nCtauErrBins));
  ws->pdf("bkgPdf")->plotOn(myPlot2_B,Name("bkgPdf"), LineColor(kBlue+2), LineWidth(2), Range("ctauErrRange"));
  ws->data("dsAB")->plotOn(myPlot2_B,Name("dataHist_Tot"), MarkerSize(.7), Binning(nCtauErrBins));
  myPlot2_B->GetYaxis()->SetRangeUser(10e-4, 10e8);
  myPlot2_B->GetXaxis()->SetTitle("l_{J/#psi} Error (mm)");
  myPlot2_B->Draw();*/

  RooDataSet* data = (RooDataSet*)ws->data("dsAB")->Clone("TMP_DATA");
  RooArgSet* cloneSet = (RooArgSet*)RooArgSet(*ws->pdf("model_A"),"model_A").snapshot(kTRUE);
  auto pdf = (RooAbsPdf*)cloneSet->find("model_A");
  pdf->setOperMode(RooAbsArg::ADirty, kTRUE);

  RooStats::SPlot sData = RooStats::SPlot("sData", "An SPlot", *data, pdf, yieldList);
  ws->import(*data, Rename("dataset_SPLOT"));
  cout<<"[INFO] Jpsi yield -> Mass Fit:"<<ws->var("nSig")->getVal()<<", sWeights :"<<sData.GetYieldFromSWeight("nSig")<<endl;
  cout<<"[INFO] Bkg  yield -> Mass Fit:"<<ws->var("nBkg")->getVal()<<", sWeights :"<<sData.GetYieldFromSWeight("nBkg")<<endl;
//create weighted data sets
//total
  TH1D* hTot = (TH1D*)ws->data("dsAB")->createHistogram(("hTot"), *ws->var("ctau3DErr"),Binning(nCtauErrBins,ctauErrLow,ctauErrHigh));
  RooDataHist* totHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hTot);
  RooHistPdf* totPdf = new RooHistPdf("totPdf","hist pdf", *ws->var("ctau3DErr"), *totHist);
//bkg
  RooDataSet* dataw_Bkg = new RooDataSet("dataw_Bkg","TMP_BKG_DATA", (RooDataSet*)ws->data("dataset_SPLOT"), 
                              RooArgSet(*ws->var("ctau3DErr"), *ws->var("nBkg_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D")), 0, "nBkg_sw");
  TH1D* hBkg = (TH1D*)dataw_Bkg->createHistogram(("hBkg"), *ws->var("ctau3DErr"),Binning(nCtauErrBins,ctauErrLow,ctauErrHigh));
  RooDataHist* bkgHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hBkg);
  RooHistPdf* bkgPdf = new RooHistPdf("bkgPdf","hist pdf", *ws->var("ctau3DErr"), *bkgHist);
//data
  RooDataSet* dataw_Sig = new RooDataSet("dataw_Sig","TMP_SIG_DATA", (RooDataSet*)ws->data("dataset_SPLOT"),
                              RooArgSet(*ws->var("ctau3DErr"), *ws->var("nSig_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D")), 0, "nSig_sw");
  TH1D* hSig = (TH1D*)dataw_Sig->createHistogram(("hSig"), *ws->var("ctau3DErr"),Binning(nCtauErrBins,ctauErrLow,ctauErrHigh));
  RooDataHist* sigHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hSig);
  RooHistPdf* sigPdf = new RooHistPdf("sigPdf","hist pdf", *ws->var("ctau3DErr"), *sigHist);
//import
  ws->import(*dataw_Sig);
  ws->import(*dataw_Bkg);
  ws->import(*totPdf);
  ws->import(*sigPdf);
  ws->import(*bkgPdf);
  cout<<ws->data("dsAB")->numEntries()<<endl;
  cout<<ws->data("dataset_SPLOT")->numEntries()<<endl;

  RooPlot* myPlot2_B = (RooPlot*)myPlot_B->Clone();
  ws->data("dataset_SPLOT")->plotOn(myPlot2_B,Name("dataHist_Tot"), MarkerSize(.7), Binning(nCtauErrBins));
  ws->pdf("totPdf")->plotOn(myPlot2_B,Name("totPdf"), LineColor(kGreen+1), Range("ctauErrRange"), LineWidth(2));
  ws->data("dataw_Sig")->plotOn(myPlot2_B,Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2), Binning(nCtauErrBins));
  ws->pdf("sigPdf")->plotOn(myPlot2_B,Name("sigPdf"),LineColor(kRed+2), LineWidth(2), Range("ctauErrRange"));
  ws->data("dataw_Bkg")->plotOn(myPlot2_B,Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue+2), MarkerColor(kBlue+2), Binning(nCtauErrBins));
  ws->pdf("bkgPdf")->plotOn(myPlot2_B,Name("bkgPdf"), LineColor(kBlue+2), LineWidth(2), Range("ctauErrRange"));
  myPlot2_B->GetYaxis()->SetRangeUser(10e-4, 10e8);
  myPlot2_B->GetXaxis()->CenterTitle();
  myPlot2_B->GetXaxis()->SetTitle("l_{J/#psi} Error (mm)");
  myPlot2_B->Draw();
  //Double_t outTot = ws->data("dsAB")->numEntries();
  //Double_t outErr = ws->data("dsAB")->reduce(Form("(ctauErr>=%.6f || ctauErr<=%.6f)", ctauErrHigh, ctauErrLow))->numEntries();
  //cout<<(outErr*100)/outTot<<endl;
  TLegend* leg_B = new TLegend(text_x+0.5,text_y-0.2,text_x+0.7,text_y); leg_B->SetTextSize(text_size);
  leg_B->SetTextFont(43);
  leg_B->SetBorderSize(0);
  leg_B->AddEntry(myPlot2_B->findObject("dataHist_Tot"),"Data","pe");
  leg_B->AddEntry(myPlot2_B->findObject("totPdf"),"Total PDF","l");
  leg_B->AddEntry(myPlot2_B->findObject("sigPdf"),"Signal","l");
  leg_B->AddEntry(myPlot2_B->findObject("bkgPdf"),"Background","l");
  leg_B->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  //drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);
  //drawText(Form("n_{J/#psi} = %.f #pm %.f",sData.GetYieldFromSWeight("nSig"),ws->var("nSig")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
  //drawText(Form("n_{Bkg} = %.f #pm %.f",sData.GetYieldFromSWeight("nBkg"), ws->var("nBkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
  cout << endl << "************** Finished SPLOT *****************" << endl << endl;
  //***********************************************************************
  //**************************** CTAU RES FIT *****************************
  //***********************************************************************
  cout << endl << "*************** Start Res Fit *****************" << endl << endl;
  RooDataSet *ctauResCutDS =(RooDataSet*)dataw_Sig->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")),*(ws->var("ctau3DErr"))),"ctau3DRes<0");
  ctauResCutDS->SetName("ctauResCutDS");
  ws->import(*ctauResCutDS);
  double entries = ws->data("ctauResCutDS")->numEntries();
  cout<<"[Info] #J/psi: "<<entries<<endl;
  ws->factory(Form("N_Jpsi[%.12f,%.12f,%.12f]", entries, entries, entries*2.0));
// create the variables for this model
  //ws->var("ctau3DRes")->setRange("ctauResRange", ctauResLow, ctauResHigh);
  int nGauss = 3;
  ws->factory("ctauRes_mean[0.0]");
  ws->factory("ctau1_CtauRes[0., -0.1, 0.1]");  ws->factory("s1_CtauRes[.5, 0., 10.]");
  ws->factory("ctau2_CtauRes[0., -0.1, 0.1]");  ws->factory("s2_CtauRes[1.12, 0., 10.]");
  ws->factory("ctau3_CtauRes[0., -0.1, 0.1]");  ws->factory("s3_CtauRes[3.37, 0., 10.]");
  ws->factory("ctau4_CtauRes[0., -0.1, 0.1]");  ws->factory("s4_CtauRes[5.37, 0., 10.]");
  ws->factory("f2_CtauRes[0.5, 0., 1.]"); ws->factory("f_CtauRes[0.5, 0., 1.]");
  ws->factory("f3_CtauRes[0.5, 0., 1.]");
// create the three PDFs
  TString varName="ctau3DRes";
  ws->factory(Form("GaussModel::%s(%s, %s, %s)", "GaussModel1_ctauRes", varName.Data(), 
                   "ctauRes_mean", //"ctau1_CtauRes",
                   "s1_CtauRes"
                   //"One", //meanSF
                   //"SF_sigma"//sigmaSF (usePerEventError?"ctauErr":"One")
                   ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s)", "GaussModel2_ctauRes", varName.Data(), 
                   "ctauRes_mean", //"ctau2_CtauRes",
                   "s2_CtauRes"
                   //"One", //meanSF
                   //"SF_sigma"//sigmaSF (usePerEventError?"ctauErr":"One")
                   ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s)", "GaussModel3_ctauRes", varName.Data(), 
                   "ctauRes_mean", //"ctau3_CtauRes",
                   "s3_CtauRes"
                   //"One", //meanSF
                   //"SF_sigma"//sigmaSF (usePerEventError?"ctauErr":"One")
                   ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s)", "GaussModel4_ctauRes", varName.Data(), 
                   "ctauRes_mean", //"ctau3_CtauRes",
                   "s4_CtauRes"
                   //"One", //meanSF
                   //"SF_sigma"//sigmaSF (usePerEventError?"ctauErr":"One")
                   ));
// combine the two PDFs
  if(nGauss==4){
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel43_ctauRes", "GaussModel4_ctauRes", "GaussModel3_ctauRes", "f3_CtauRes"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel23_ctauRes", "GaussModel43_ctauRes", "GaussModel2_ctauRes", "f2_CtauRes"));}
  else if(nGauss==3){
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel23_ctauRes", "GaussModel3_ctauRes", "GaussModel2_ctauRes", "f2_CtauRes"));}
  ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModelCOND_ctauRes", "GaussModel1_ctauRes", "GaussModel23_ctauRes", "f_CtauRes"));
  //GaussModelCOND_ctauRes is total model with 3 gm..
  //ws->import(GaussModelCOND_ctauRes); //wrong
//#1
  //RooProdPdf pdf1("GaussModel_ctauRes", "", *ws->pdf("sigPdf"), Conditional(*ws->pdf("GaussModelCOND_ctauRes"), RooArgList(*ws->var("N_Jpsi"))));
  //ws->import(pdf1);
  //ws->factory("RooExtendPdf::GaussModel_Tot_ctauRes(GaussModel_ctauRes,N_Jpsi)");
//#2
  //RooAddPdf* ctauResModel = new RooAddPdf();
  //ctauResModel = new RooAddPdf("GaussModel_Tot","Jpsi Res",RooArgList("GaussModelCOND_ctauRes"),RooArgList(*nSig));
//#3
  RooAbsPdf *ctauResModel = ctauResModel = new RooAddPdf("GaussModel_Tot", "GaussModel_Tot", *ws->pdf("GaussModelCOND_ctauRes"), *ws->var("N_Jpsi"));
  ws->import(*ctauResModel);
  //setFixedVarsToContantVars(ws);
  //void setFixedVarsToContantVars(RooWorkspace& ws)
  //{
  //RooArgSet listVar = ws->allVars();
  //TIterator* parIt = listVar.createIterator();
  //for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
  //  if ( it->getMin()==it->getMax() && !it->isConstant() ) {
  //    cout << "[INFO] Setting " << it->GetName() << " constant!" << endl;
  //    it->setConstant(kTRUE);
  //  }
  //}
  //};
  c_C->cd();
  c_C->SetLogy();
  RooPlot* myPlot2_C = (RooPlot*)myPlot_C->Clone();
  bool isWeighted = ws->data("ctauResCutDS")->isWeighted();
  RooFitResult* fitCtauRes = ws->pdf("GaussModel_Tot")->fitTo(*ctauResCutDS, Save(), Range(ctauResLow,0), SumW2Error(isWeighted), Extended(kTRUE), NumCPU(12), PrintLevel(3));
  ws->data("ctauResCutDS")->plotOn(myPlot2_C, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0), 
      MarkerSize(.7), Binning(120), LineColor(kRed+2), MarkerColor(kRed+2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_Tot"), Precision(1e-5), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), LineColor(kBlack));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm1"), Precision(1e-5), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel1_ctauRes")), LineColor(kGreen+2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm2"), Precision(1e-5), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel2_ctauRes")), LineColor(kRed+2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm3"), Precision(1e-5), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel3_ctauRes")), LineColor(kBlue+2));
  if(nGauss==4){ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm4"), Precision(1e-5), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel4_ctauRes")), LineColor(kMagenta+2));}
  myPlot2_C->GetYaxis()->SetRangeUser(10e-4, 10e6);
  myPlot2_C->GetXaxis()->CenterTitle();
  myPlot2_C->GetXaxis()->SetTitle("#frac{l_{J/#psi}}{#sigma_{J/#psi}}");
  myPlot2_C->Draw();
  TLegend* leg_C = new TLegend(text_x+0.5,text_y-0.27,text_x+0.65,text_y-0.03); leg_C->SetTextSize(text_size);
  leg_C->SetTextFont(43);
  leg_C->SetBorderSize(0);
  leg_C->AddEntry(myPlot2_C->findObject("dataHist_ctauRes"),"Data_Sig","pe");
  leg_C->AddEntry(myPlot2_C->findObject("modelHist_Tot"),"Total PDF","l");
  leg_C->AddEntry(myPlot2_C->findObject("modelHist_gm1"),"Gauss 1","l");
  leg_C->AddEntry(myPlot2_C->findObject("modelHist_gm2"),"Gauss 2","l");
  leg_C->AddEntry(myPlot2_C->findObject("modelHist_gm3"),"Gauss 3","l");
  if(nGauss==4)leg_C->AddEntry(myPlot2_C->findObject("modelHist_gm4"),"Gauss 4","l");
  leg_C->Draw("same");
  cout<<"s2/s1: "<<ws->var("s2_CtauRes")->getVal()/ws->var("s1_CtauRes")->getVal()<<endl;
  cout<<"s3/s2: "<<ws->var("s3_CtauRes")->getVal()/ws->var("s2_CtauRes")->getVal()<<endl;
  cout<<"s1: "<<ws->var("s1_CtauRes")->getVal()<<endl;
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);
  
  //drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError()),text_x+0.5,text_y,text_color,text_size);
  drawText(Form("N_{J/#psi} = %.f", entries ),text_x+0.5,text_y,text_color,text_size);
  //drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c", ws->ws->var("s2_CtauRes")->getVal()/ws->var("s1_CtauRes")->getVal(), ptHigh ),text_x+0.5,text_y,text_color,text_size);
  cout << endl << "************* Finished Sig Res Fit ****************" << endl << endl;
  cout << endl << "************ Start MC Ctau True Fit ***************" << endl << endl;
//MC NP ctau true
  double entries_True = wsmc->data("reducedDS_MC")->numEntries();
  wsmc->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, entries_True, entries*1.1));
  wsmc->factory("lambdaDSS[0.3, 0., 1.]");
  // create the PDF
  wsmc->factory(Form("TruthModel::%s(%s)", "TruthModel_ctauTrue", "ctau3Dtrue"));
  wsmc->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "ctauTruePdf", "ctau3Dtrue",  
        "lambdaDSS", 
        "TruthModel_ctauTrue"));

  RooAbsPdf *ctauTrueModel = ctauTrueModel = new RooAddPdf("TrueModel_Tot", "TrueModel_Tot", *wsmc->pdf("ctauTruePdf"), *wsmc->var("N_Jpsi_MC"));
  wsmc->import(*ctauTrueModel);
  
  //RooAbsPdf *themodel = new RooAddPdf("pdfNoPR", "pdfNoPR", pdfList);
  //wsmc->import(*themodel);
  //wsmc->pdf(pdfName.c_str())->setNormRange("CtauTrueWindow");
  c_D->cd();
  c_D->SetLogy();
  RooFitResult* fitCtauTrue = wsmc->pdf("TrueModel_Tot")->fitTo(*reducedDS_MC, Save(), Range("ctauTrueRange"), Extended(kTRUE), NumCPU(12), PrintLevel(-1));
  RooPlot* myPlot2_D = (RooPlot*)myPlot_D->Clone();
  wsmc->data("reducedDS_MC")->plotOn(myPlot2_D,Name("MCHist_Tot"), MarkerSize(.7), Binning(nCtauBins));
  wsmc->pdf("TrueModel_Tot")->plotOn(myPlot2_D,Name("MCpdf_Tot"), LineColor(kRed+2));
  myPlot2_D->GetYaxis()->SetRangeUser(10e-2, wsmc->data("reducedDS_MC")->sumEntries()*10);
  myPlot2_D->GetXaxis()->SetRangeUser(-1, 7);
  myPlot2_D->GetXaxis()->CenterTitle();
  myPlot2_D->GetXaxis()->SetTitle("l_{J/#psi} MC True (mm)");
  myPlot2_D->Draw();
  
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),text_x+0.5,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x+0.5,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x+0.5,text_y-y_diff,text_color,text_size);
  //***********************************************************************
  //**************************** Bkg CTAU FIT *****************************
  //***********************************************************************
  cout << endl << "************** Start Ctau Fit *****************" << endl << endl;
//
  ws->factory("zeroMean[0.]");
  ws->factory("SF_sigma[0.5, 0., 1.0]");
  ws->factory("lambdaDSS_Bkg[0.03, 0., 1.]");
  ws->factory("lambdaDF[0.01, 0., 0.1]");
  ws->factory("lambdaDDS[0.01, 0., 0.1]");
  ws->factory("fDFSS[0.5, 0., 1.]");
  ws->factory("fDLIV[0.5, 0., 1.]");

  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes1", "ctau3D", 
                   "ctauRes_mean", //"ctau1_CtauRes",
                   "s1_CtauRes",
                   "zeroMean",
                   "ctau3DErr"
                   ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes2", "ctau3D", 
                   "ctauRes_mean", //"ctau2_CtauRes",
                   "s2_CtauRes",
                   "zeroMean",
                   "ctau3DErr"
                   ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes3", "ctau3D", 
                   "ctauRes_mean", //"ctau3_CtauRes",
                   "s3_CtauRes",
                   "zeroMean",
                   "ctau3DErr"
                   ));

  ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "ctauRes32", "ctauRes3", "ctauRes2", "f2_CtauRes"));
  ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "ctauResCond", "ctauRes32", "ctauRes1", "f_CtauRes"));

  //RooArgSet* cloneSet1 = (RooArgSet*)RooArgSet(*ws->pdf("GaussModelCOND_ctauRes"),"GaussModelCOND_ctauRes").snapshot(kTRUE);
  //auto ctauResCond = (RooAbsPdf*)cloneSet1->find("GaussModelCOND_ctauRes");
  //ctauResCond->setOperMode(RooAbsArg::ADirty, kTRUE);

  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS", "ctau3D",
                 "lambdaDSS_Bkg", "ctauResCond"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF", "ctau3D",
                 "lambdaDF", "ctauResCond"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", "pdfCTAUDDS", "ctau3D",
                  "lambdaDDS", "ctauResCond"));
  
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU1",
                  "fDFSS", "pdfCTAUDSS", "pdfCTAUDF"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND",
                  "fDLIV", "pdfCTAU1", "pdfCTAUDDS"));

  RooAbsPdf *ctauBkgModel = ctauBkgModel = new RooAddPdf("BkgModel_Tot", "BkgModel_Tot", *ws->pdf("pdfCTAUCOND"), *ws->var("nBkg"));
  ws->import(*ctauBkgModel);

  c_E->cd();
  c_E->SetLogy();
  RooPlot* myPlot2_E = (RooPlot*)myPlot_E->Clone();
  RooFitResult* fitCtauBkg = ws->pdf("BkgModel_Tot")->fitTo(*dataw_Bkg, Save(), Range("ctauRange"), Extended(kTRUE), NumCPU(12), PrintLevel(3));
  //ws->data("dataw_Sig")->plotOn(myPlot2_E,Name("dataHist_Tot"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kRed+2), MarkerColor(kRed+2));
  ws->data("dataw_Bkg")->plotOn(myPlot2_E,Name("dataHist_ctauBkg"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kBlue+2), MarkerColor(kBlue+2));
  //ws->pdf("BkgModel_Tot")->plotOn(myPlot2_E,Name("ctauBkg_Tot"), LineColor(kBlack));
  ws->pdf("BkgModel_Tot")->plotOn(myPlot2_E,Name("BKG"), 
      Normalization(ws->data("dataw_Bkg")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE),
      FillStyle(1001), FillColor(kAzure-9), LineColor(kBlack), LineWidth(3), DrawOption("LCF"), Precision(1e-4));
      //FillStyle(1001), FillColor(kAzure-9), LineColor(kBlack), VLines(), DrawOption("LCF"), Precision(1e-4));
  ws->pdf("BkgModel_Tot")->plotOn(myPlot2_E,Name("BKG"), 
      Normalization(ws->data("dataw_Bkg")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE), Components("ctauResCond"),
      LineWidth(2), LineColor(kRed+2));
  //ws->pdf("BkgModel_Tot")->plotOn(myPlot2_E, Name("test"), Components(*ws->pdf("pdfCTAU1")), LineColor(kRed+2));
  //ws->pdf("BkgModel_Tot")->plotOn(myPlot2_E,Name("resPdf"), Components(*ws->pdf("pdfCTAUDSS")), LineColor(kRed+2));
  //ws->pdf("BkgModel_Tot")->plotOn(myPlot2_E,Name("resPdf"), Components(*ws->pdf("pdfCTAUDF")), LineColor(kBlue+2));
  //ws->pdf("BkgModel_Tot")->plotOn(myPlot2_E,Name("resPdf"), Components(*ws->pdf("pdfCTAUDDS")), LineColor(kMagenta+2));
  ws->data("dataw_Bkg")->plotOn(myPlot2_E,Name("dataHist_ctauBkg"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kBlue+2), MarkerColor(kBlue+2));
  myPlot2_E->GetYaxis()->SetRangeUser(10e-2, 10e6);
  myPlot2_E->GetXaxis()->SetRangeUser(-4, 6);
  myPlot2_E->GetXaxis()->SetTitle("l_{J/#psi} (mm)");
  myPlot2_E->Draw();
  TLegend* leg_E = new TLegend(text_x+0.5,text_y-0.17,text_x+0.65,text_y-0.03); leg_E->SetTextSize(text_size);
  leg_E->SetTextFont(43);
  leg_E->SetBorderSize(0);
  leg_E->AddEntry(myPlot2_E->findObject("dataHist_ctauBkg"),"Data_Bkg","pe");
  leg_E->AddEntry(myPlot2_E->findObject("ctauBkg_Tot"),"Total PDF","fl");
  //leg_E->AddEntry(myPlot2_E->findObject("test"),"?? PDF","l");
  leg_E->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);

//  if ( ws->pdf("pdfCTAUMASS_JpsiNoPR") ) {
//    ws->pdf(pdfTotName.c_str())->plotOn(frame,Name("JPSINOPR"),
//        Components(RooArgSet(*ws->pdf("pdfCTAUMASS_JpsiNoPR"), *ws->pdf("pdfCTAUMASS_Bkg"))),
//        ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsAB"), kTRUE),
//        Normalization(dsAB, RooAbsReal::NumEvent),
//        LineColor(kGreen+3), LineStyle(1), Precision(1e-4), NumCPU(32)
//        );
//  }

  //***********************************************************************
  //*************************** DRAW CTAU FIT *****************************
  //***********************************************************************
  c_F->cd();
  c_F->SetLogy();
  double normBkg = ws->data("dsAB")->sumEntries()/ws->data("dataw_Bkg")->sumEntries();
  cout<<"[Info]: "<<normBkg<<endl;
  RooPlot* myPlot2_F = (RooPlot*)myPlot_F->Clone();
  myPlot2_F->updateNormVars(RooArgSet(*ws->var("mass"), *ws->var("ctau3D"), *ws->var("ctau3DErr")));
  ws->pdf("BkgModel_Tot")->plotOn(myPlot2_F,Name("BKG"), 
      //Normalization(normBkg, RooAbsReal::NumEvent), NormRange("ctauRange"),
      Normalization(ws->data("dataw_Bkg")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE),
      FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("CF"), Precision(1e-4));
  ws->data("dsAB")->plotOn(myPlot2_F,Name("dataHist_ctauTot"), DataError(RooAbsData::SumW2), XErrorSize(0));
  //ws->pdf("BkgModel_Tot")->plotOn(myPlot2_F,Name("resPdf"), Components(*ws->pdf("pdfCTAUDSS")), LineColor(kRed+2));
  //ws->pdf("BkgModel_Tot")->plotOn(myPlot2_F,Name("resPdf"), Components(*ws->pdf("pdfCTAUDF")), LineColor(kBlue+2));
  //ws->pdf("BkgModel_Tot")->plotOn(myPlot2_F,Name("resPdf"), Components(*ws->pdf("pdfCTAUDDS")), LineColor(kMagenta+2));
  myPlot2_F->GetYaxis()->SetRangeUser(10e-2, 10e6);
  myPlot2_F->GetXaxis()->SetRangeUser(-4, 6);
  myPlot2_F->GetXaxis()->SetTitle("l_{J/#psi} (mm)");
  myPlot2_F->Draw();
  TLegend* leg_F = new TLegend(text_x+0.5,text_y-0.17,text_x+0.65,text_y-0.03); leg_F->SetTextSize(text_size);
  leg_F->SetTextFont(43);
  leg_F->SetBorderSize(0);
  leg_F->AddEntry(myPlot2_F->findObject("dataHist_ctauTot"),"Data","pe");
  leg_F->AddEntry(myPlot2_F->findObject("BKG"),"Background PDF","fl");
  //leg_F->AddEntry(myPlot2_F->findObject("test"),"?? PDF","l");
  leg_F->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);

  c_A->Update();
  c_B->Update();
  c_C->Update();
  c_D->Update();
  c_E->Update();
  c_F->Update();

  TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);
  TString outFileName;
  if (whichModel){
    outFileName = Form("test_Sim_BSplit_altfitresults_Jpsi_%s.root",kineLabel.Data());
  }
  else {
    outFileName = Form("test_Sim_BSplit_nomfitresults_Jpsi_%s.root",kineLabel.Data());
  }
  TFile* outf = new TFile(outFileName,"recreate");
  ws->Write();
  outf->Close();

} 
