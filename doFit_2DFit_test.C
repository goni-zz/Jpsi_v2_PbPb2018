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
       float ptLow=6.5, float ptHigh=50, 
       float yLow=0, float yHigh=2.4,
       int cLow=0, int cHigh=200,
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
  
  double ctauErrLow = 0., ctauErrHigh = 0.18;
  double binWidth = 0.0025;

  TFile* f1; TFile* f2;
  TString kineCut;
  TString SigCut;
  TString BkgCut;

  //Select Data Set
    //f1 = new TFile("skimmedFiles/OniaRooDataSet_isMC0_JPsi.root");
    f1 = new TFile("skimmedFiles/OniaRooDataSet_isMC0_JPsi1SW_2020819.root");
    f2 = new TFile("skimmedFiles/OniaRooDataSet_isMC1_JPsi1SW_2020821.root");
    //kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.5",ptLow, ptHigh, yLow, yHigh);
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
  RooWorkspace *wsmc = new RooWorkspace("workspaceMC");
  ws->import(*dataset);
  wsmc->import(*datasetMC);
  cout << "####################################" << endl;
  RooDataSet *reducedDS_A = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), SigCut.Data() );
  RooDataSet *reducedDS_B = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), BkgCut.Data() );
  RooDataSet *reducedDS_MC = (RooDataSet*)datasetMC->reduce(RooArgSet(*(wsmc->var("ctau3Dtrue")), *(wsmc->var("mass")), *(wsmc->var("pt")), *(wsmc->var("y"))), kineCut.Data() );
  reducedDS_A->SetName("reducedDS_A");
  reducedDS_B->SetName("reducedDS_B");
  reducedDS_MC->SetName("reducedDS_MC");
  //ws->import(*reducedDS_A);
  //ws->import(*reducedDS_B);
  wsmc->import(*reducedDS_MC);
  //ws->var("mass")->setRange(massLow, massHigh);
  //ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  //ws->var("ctau3DErr")->setRange(ctauLow, ctauHigh);
  //ws->var("ctau3DRes")->setRange(ctauResLow, ctauResHigh);
  //reducedDS_A->Print();
  //reducedDS_B->Print();
  reducedDS_MC->Print();
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

  //***********************************************************************
  //****************************** MASS FIT *******************************
  //***********************************************************************
  //         The order is {sigma_1,  x, alpha_1, n_1,   f, err_mu, err_sigma, m_lambda}
  double paramsupper[8] = {0.2,    3.0,   3.321, 5.0, 1.0,   15.0,      15.0,     25.0};
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
  RooRealVar *nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",10000,0,5000000); 

  //Build the model
//Model A: Mass
  RooAddPdf* model_A = new RooAddPdf();
  model_A = new RooAddPdf("model_A","Jpsi + Bkg",RooArgList(*cb_A, *bkg),RooArgList(*nSig,*nBkg));
  //model_A = new RooAddPdf("model_A","PR Jpsi + NP Jpsi + Bkg",RooArgList(*cb_1_A, *cb_2_A, *bkg),RooArgList(*nSigPR,*nSigNP,*nBkg));
  ws->import(*model_A);
  cout << endl << "********* Starting Mass Dist. Fit **************" << endl << endl;
  RooFitResult* fitResSim = ws->pdf("model_A")->fitTo(*dsAB,Save(), Hesse(kTRUE), Range(massLow,massHigh), Timer(kTRUE), Extended(kTRUE), NumCPU(12));
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
  myPlot2_A->GetYaxis()->SetLabelSize(0.054);
  //myPlot2_A->GetYaxis()->SetRangeUser(2*10, 2*10e3);
  myPlot2_A->GetYaxis()->SetRangeUser(ws->var("nSig")->getVal()/100, ws->var("nSig")->getVal());
  //myPlot2_A->SetMinimum(2*10);
  //myPlot2_A->GetXaxis()->SetLabelSize(0);
  myPlot2_A->GetXaxis()->SetRangeUser(massLow,massHigh);
  //myPlot2_A->GetXaxis()->SetTitleSize(0);
  myPlot2_A->Draw();

  float text_x = 0.15;
  float text_y = 0.816;
  float y_diff = 0.05;
  float text_size = 19;
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
  myPlot2_B->Draw();
  //Double_t outTot = ws->data("dsAB")->numEntries();
  //Double_t outErr = ws->data("dsAB")->reduce(Form("(ctauErr>=%.6f || ctauErr<=%.6f)", ctauErrHigh, ctauErrLow))->numEntries();
  //cout<<(outErr*100)/outTot<<endl;
  //drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  //if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  //else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  //drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);
  //drawText(Form("n_{J/#psi} = %.f #pm %.f",sData.GetYieldFromSWeight("nSig"),ws->var("nSig")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
  //drawText(Form("n_{Bkg} = %.f #pm %.f",sData.GetYieldFromSWeight("nBkg"), ws->var("nBkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
  cout << endl << "************** Finished SPLOT *****************" << endl << endl;
  //***********************************************************************
  //**************************** CTAU RES FIT *****************************
  //***********************************************************************
  cout << endl << "*************** Start Res Fit *****************" << endl << endl;
// create the variables for this model
  ws->var("ctau3DRes")->setRange("ctauResRange", ctauResLow, 0);
  ws->factory("One[1.0]");
  ws->factory("ctau1_CtauRes[0,-1.,1.]");  ws->factory("s1_CtauRes[0.012, 0., 1.]");
  ws->factory("ctau2_CtauRes[0,-1.,1.]");  ws->factory("s2_CtauRes[0.012, 0., 1.]");
  ws->factory("ctau3_CtauRes[0,-1.,1.]");  ws->factory("s3_CtauRes[0.012, 0., 1.]");
  ws->factory("f2_CtauRes[0.5, 0., 1.]");ws->factory("f_CtauRes[0.5, 0., 1.]");
// create the three PDFs
  TString varName="ctau3DRes";
  ws->factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", "GaussModel1_ctauRes", varName.Data(), "ctau1_CtauRes", "s1_CtauRes",
                   "ctau3DErr"//sigmaSF (usePerEventError?"ctauErr":"One")
                   ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", "GaussModel2_ctauRes", varName.Data(), "ctau2_CtauRes", "s2_CtauRes",
                   "ctau3DErr"//sigmaSF (usePerEventError?"ctauErr":"One")
                   ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", "GaussModel3_ctauRes", varName.Data(), "ctau3_CtauRes", "s3_CtauRes",
                   "ctau3DErr"//sigmaSF (usePerEventError?"ctauErr":"One")
                   ));
// combine the two PDFs
  ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel23_ctauRes", "GaussModel2_ctauRes", "GaussModel3_ctauRes", "f2_CtauRes"));
  ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModelCOND_ctauRes", "GaussModel1_ctauRes", "GaussModel23_ctauRes", "f_CtauRes"));
  //ws->import(GaussModelCOND_ctauRes);
  //RooProdPdf pdf1("GaussModel_ctauRes", "", *ws->pdf("sigPdf"), Conditional(*ws->pdf("GaussModelCOND_ctauRes"), RooArgList(*ws->var("ctau3DRes"))));
  //ws->import(pdf1);
  //ws->factory("RooExtendPdf::GaussModelTot_ctauRes(GaussModel_ctauRes,nSig_sw)");
  RooAbsPdf *themodel = themodel = new RooAddPdf("GaussModel_Tot", "GaussModel_Tot", *ws->pdf("GaussModelCOND_ctauRes"), *ws->var("nSig_sw"));
  ws->import(*themodel);
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
  //bool isWeighted = ws->data(*dsName2FitCut.c_str())->isWeighted();
  RooFitResult* fitCtauRes = ws->pdf("GaussModel_Tot")->fitTo(*dataw_Sig, Save(), Extended(kTRUE), NumCPU(12));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_Tot"), Normalization(ws->data("dataw_Sig")->sumEntries(), RooAbsReal::NumEvent));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm1"), Components(*ws->pdf("GaussModel1_ctauRes")), Normalization(ws->data("dataw_Sig")->sumEntries(), RooAbsReal::NumEvent), LineColor(kGreen+2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm2"), Components(*ws->pdf("GaussModel2_ctauRes")), Normalization(ws->data("dataw_Sig")->sumEntries(), RooAbsReal::NumEvent), LineColor(kRed+2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm3"), Components(*ws->pdf("GaussModel3_ctauRes")), Normalization(ws->data("dataw_Sig")->sumEntries(), RooAbsReal::NumEvent), LineColor(kBlue+2));
  ws->data("dataw_Sig")->plotOn(myPlot2_C,Name("dataHist_Sig"), MarkerSize(.7), Binning(nCtauBins), CutRange("ctauResRange"), DataError(RooAbsData::SumW2), XErrorSize(0), LineColor(kRed+2), MarkerColor(kRed+2));
  ws->data("dataw_Bkg")->plotOn(myPlot2_C,Name("dataHist_Bkg"), MarkerSize(.7), Binning(nCtauBins), CutRange("ctauResRange"), DataError(RooAbsData::SumW2), XErrorSize(0), LineColor(kBlue+2), MarkerColor(kBlue+2));
  myPlot2_C->GetYaxis()->SetRangeUser(10e-4, 10e6);
  myPlot2_C->Draw();

  cout << endl << "************* Finished Sig Res Fit ****************" << endl << endl;

  cout << endl << "************** Start Ctau Fit *****************" << endl << endl;
//MC NP ctau true
  c_D->cd();
  c_D->SetLogy();
  RooPlot* myPlot2_D = (RooPlot*)myPlot_D->Clone();
  wsmc->data("reducedDS_MC")->plotOn(myPlot2_D,Name("MCHist_Tot"), MarkerSize(.7), Binning(nCtauBins));

  cout << endl << "************** Start Ctau Fit *****************" << endl << endl;
  
  c_E->cd();
  c_E->SetLogy();
  RooPlot* myPlot2_E = (RooPlot*)myPlot_E->Clone();
  ws->data("dataw_Sig")->plotOn(myPlot2_E,Name("dataHist_Tot"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kRed+2), MarkerColor(kRed+2));
  ////ws->pdf("model_ct")->plotOn(myPlot2_D,Name("totPdf"), LineColor(kGreen+1), Range("ctauRange"), LineWidth(2));
  //ws->data("reducedDS_A")->plotOn(myPlot2_D,Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2), Binning(nCtauBins));
  ////ws->pdf("sigPdf")->plotOn(myPlot2_D,Name("sigPdf"),LineColor(kRed+2), LineWidth(2), Range("ctauRange"));
  //ws->data("reducedDS_B")->plotOn(myPlot2_D,Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue), MarkerColor(kBlue+2), Binning(nCtauBins));
  ////ws->pdf("bkgPdf")->plotOn(myPlot2_D,Name("bkgPdf"), LineColor(kBlue+2), LineWidth(2), Range("ctauRange"));
  ////ws->data("dsAB")->plotOn(myPlot2_D,Name("dataHist_Tot"), MarkerSize(.7), Binning(nCtauBins));
  //myPlot2_D->GetYaxis()->SetRangeUser(10e-4, 10e6);
  ////myPlot2_D->GetXaxis()->SetTitle("l_{J/#psi} (mm)");
  //myPlot2_D->Draw();
//B:Ctau distribution with 3 decay pdf
  //RooRealVar tau("tau","tau", 1.548); // some good start value
  //RooTruthModel idealres("idealres","Ideal resolution model", *ws->var("ctau3D"));//needed MC?
  //RooDecay decay("decay", "decay", *ws->var("ctau3D"), tau, idealres, RooDecay::SingleSided);

  //RooRealVar dt("dt","per-event error on t",0.01,5) ;
  //RooRealVar bias("bias","bias",0,-3,5) ;
  //RooRealVar sigma("sigma","per-event error scale factor",1,0.1,2) ; 
  //RooGaussModel gm("gm1","gauss model scaled by per-event error",*ws->var("ctau3D"),bias,sigma);

  //RooRealVar testsig("testsig","testsig",ws->var("nSig")->getVal(),ws->var("nSig")->getVal(),ws->var("nSig")->getVal());
  //RooAddPdf* ctau_sig = new RooAddPdf("ctau","ctau",RooArgSet(decay),RooArgList(testsig));
//B: Ctau Bkg function
  //RooRealVar l_DDS("l_DDS","l_DDS", 0.1,  0., 1.); 
  //RooRealVar l_DF("l_DF","l_DF", 0.1,  0., 1.); 
  //RooRealVar l_DSS("l_DSS","l_DSS", 0.1,  0., 1.); 
  //
  //RooRealVar fDLIV("fDLIV","fDLIV",0.1, 0, 1.) ;
  //RooRealVar fDFSS("fDFSS","fDFSS",0.1, 0, 1.) ;
  //RooRealVar bbkg( "bbkg", "bbkg" ,0.5, 0, 1.) ;
  //RooGenericPdf *DSS = new RooGenericPdf("DSS","Background",
      //"@6*(@4*(@5*TMath::Exp(-abs(@1)*@0)+(@5*TMath::Exp(abs(@2)*@0)))+(1-@4)*TMath::Exp(-abs(@3)*@0))+(1-@6)*@7",RooArgList(*(ws->var("ctau3D")),l_DDS,l_DF,l_DSS,fDLIV,fDFSS,bbkg,*(ws->var("ctau3DErr"))));
      //"(@4*(@5*TMath::Exp(-abs(@1)*@0)+((1-@5)*TMath::Exp(abs(@2)*@0)))+(1-@4)*TMath::Exp(-abs(@3)*@0))",RooArgList(*(ws->var("ctau3D")),l_DDS,l_DF,l_DSS,fDLIV,fDFSS));
  //0:ctau, 1:l_DDS, 2:l_DF, 3:l_DSS, 4:fDLIV, 5:fDFSS
  
  //RooRealVar ctmean("ctmean","Mean of Gaussian",0.0,-0.3,0.3);
  //RooGaussian ctgauss("ctgauss","gauss(x,mean,sigma)",*ws->var("ctau3D"),ctmean,*ws->var("ctau3DErr")) ;

  //RooFormulaVar ff("ff", "bbkg*(fDLIV*(fDFSS*DSS)+(1-fDLIV)*DF)+(1-fDLIV)*DDS+(1-bbkg)*ctgauss", RooArgList(*ws->var("ctau3D"),fDLIV,fDLIV,bbkg,ctgauss, *DSS, *DF, *DDS));
  //RooFormulaVar ff1("ff1", "bbkg*(fDLIV*(fDFSS*DSS)", RooArgList(fDLIV,fDFSS,bbkg,*DSS));
  //RooFormulaVar ff2("ff2", "bbkg*((1-fDLIV)*DF)", RooArgList(fDLIV,bbkg,*DF));
  //RooFormulaVar ff3("ff3", "(1-fDLIV)*DDS", RooArgList(fDLIV,*DDS));
  //RooFormulaVar ff4("ff4", "(1-bbkg)*ctgauss", RooArgList(bbkg,ctgauss));
  //RooFormulaVar ff("ff", "@2*(@0*(@1*@6)+(1-@1)*@5)+(1-@0)*@4+(1-@2)*@3", RooArgList(fDLIV,fDFSS,bbkg,ctgauss, *DSS, *DF, *DDS));
  //RooAddPdf* ctau_bkg = new RooAddPdf("ctau_bkg","ctau_bkg",RooArgList(*DSS,ctgauss), RooArgList(bbkg));/////used
  //RooAddPdf* ctau_bkg = new RooAddPdf("ctau_bkg","ctau_bkg",RooArgList(ctgauss, *DSS, *DF, *DDS),RooArgList(*ws->var("ctau3D"),fDLIV,fDLIV,bbkg));
  //RooAddPdf* ctau_bkg = new RooAddPdf("ctau_bkg","ctau_bkg",RooArgList(ctgauss),RooArgList(*ws->var("ctau3D")));
  //RooAddPdf::ctau_bkg[ bbkg * (fDLIV*(fDFSS*DSS)+(1-fDLIV)*DF)+(1-fDLIV)*DDS+(1-bbkg)* (*ws->var("ctau3DErr")) ];

  //RooFormulaVar fracAlpha("fracAlpha","@0/(@0+@1)",RooArgList(*nBkg,*nSig)); ws->import(fracAlpha);

  //RooRealVar testbkg("testbkg","testbkg",ws->var("nBkg")->getVal(),ws->var("nBkg")->getVal()*0.9,ws->var("nBkg")->getVal());
//Model B: Ctau
  //RooAddPdf* model_B = new RooAddPdf();
  //model_B = new RooAddPdf("model_B","Bkg",RooArgList(*ctau_bkg), RooArgList(testbkg));
  //model_B = new RooAddPdf("model_B","Jpsi + Bkg",RooArgList(*ctau_sig, *ctau_bkg), RooArgList(testsig, testbkg));
  //model_B = new RooAddPdf("model_B","Jpsi + Bkg",RooArgList(*ctPRRes, *CkBkgTot),RooArgList(*nSig, *nBkg));
  //ws->import(*model_B);
  // Construct simultaneous PDF
  //RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simPdf", tp);
  //simPdf->addPdf(*(ws->pdf("model_B")),"A") ;
  //simPdf->addPdf(*(ws->pdf("model_B")),"B") ;
  //ws->import(*simPdf);

  //RooFitResult* fitResSim_ctau = ws->pdf("model_B")->fitTo(*reducedDS_B ,Save(), Hesse(kTRUE), Range(ctauLow, ctauHigh), Timer(kTRUE));
  //RooFitResult* fitResSim_ctau = ws->pdf("model_B")->fitTo(*dsAB ,Save(), Hesse(kTRUE), Range(ctauLow, ctauHigh), Timer(kTRUE), Extended(kTRUE));
  //ws->pdf("model_B")->plotOn(myPlot2_B,Name("modelHist_B"), LineColor(kBlack));
  //ws->pdf("model_B")->plotOn(myPlot2_B,Name("Sig_B"),Components(RooArgSet(*ctau_sig)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  //ws->pdf("model_B")->plotOn(myPlot2_B,Name("bkgPDF_B"),Components(RooArgSet(*ctau_bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  //myPlot2_B->SetFillStyle(4000);
  //myPlot2_B->GetYaxis()->SetTitleOffset(1.43);
  //myPlot2_B->GetYaxis()->CenterTitle();
  //myPlot2_B->GetYaxis()->SetTitleSize(0.058);
  //myPlot2_B->GetYaxis()->SetLabelSize(0.054);
  //myPlot2_B->GetXaxis()->SetLabelSize(0);
  //myPlot2_B->GetXaxis()->SetRangeUser(ctauLow,ctauHigh);
  //myPlot2_B->GetYaxis()->SetRangeUser(10e-2, 10e5);
  //myPlot2_B->GetXaxis()->SetTitleSize(0);
  //myPlot2_B->Draw();

  c_A->Update();
  c_B->Update();
  c_C->Update();
  c_D->Update();
  c_E->Update();

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
