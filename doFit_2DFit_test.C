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

void pullDist(RooWorkspace* ws, TPad* Pad, TCanvas* canvas, RooPlot* frame, RooHist* hpull, string dataHist, string modelHist, string variable, int nBins, float Low, float High, string titleX);
void printChi2(RooWorkspace* myws, TPad* Pad, RooPlot* frame, RooFitResult* fitRes, string varLabel, string dataLabel, string pdfLabel, int nBins, bool useDefaultName=true);
void printChi2_test(RooWorkspace* myws, TPad* Pad, RooPlot* frame, RooHist* hpull, RooFitResult* fitRes, string varLabel, string dataLabel, string pdfLabel, int nBins, bool useDefaultName=true);

void doFit_2DFit_test( 
       float ptLow=6.5, float ptHigh=50, 
       float yLow=0, float yHigh=0.6,
       int cLow=0, int cHigh=200,
       float muPtCut=0.0,
       bool whichModel=0,   // Nominal = 0. Alternative = 1.
       int ICset = 1
			) 
{
  gStyle->SetEndErrorSize(0);

  float massLow = 2.6, massHigh = 3.5;
  int   nMassBin  = 36; //(massHigh-massLow)*30;

  float ctauLow = -4, ctauHigh = 6.5;
  float ctauResLow = -20, ctauResHigh = 20;
  int   nCtauErrBins = 72;
  int   nCtauResBins = 72;
  int   nCtauBins  = (ctauHigh-ctauLow)*10;
  int   nCtauTrueBins  = 72;
  
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
    //TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

  kineCut = kineCut;
  SigCut =  SigCut;
  BkgCut =  BkgCut;
  
  //kineCut = kineCut;
  //SigCut = SigCut;
  //BkgCut = BkgCut;
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
  ws->import(*reducedDS_MC);
  reducedDS_MC->Print();
  ws->var("ctau3Dtrue")->setRange(0, 5);
  ws->var("ctau3Dtrue")->setRange("ctauTrueRange", -1, 6.5);
  ws->var("ctau3Dtrue")->Print();
  
  TCanvas* c_A =  new TCanvas("canvas_A","My plots",4,4,550,520);
  c_A->cd();
  TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0, 0.16, 0.98, 1.0);
  pad_A_1->SetTicks(1,1);
  pad_A_1->Draw(); pad_A_1->cd();
  RooPlot* myPlot_A = ws->var("mass")->frame(nMassBin); // bins
  myPlot_A->SetTitle("");
  ws->data("dsAB")->plotOn(myPlot_A,Name("dataHist_A"));

  TCanvas* c_B =  new TCanvas("canvas_B","My plots",554,4,550,520);
  c_B->cd();
  TPad *pad_B_1 = new TPad("pad_B_1", "pad_B_1", 0, 0.16, 0.98, 1.0);
  pad_B_1->SetTicks(1,1);
  pad_B_1->Draw(); pad_B_1->cd();
  RooPlot* myPlot_B = ws->var("ctau3DErr")->frame(nCtauErrBins); // bins
  myPlot_B->SetTitle("");
 // ws->data("dsAB")->plotOn(myPlot_B,Name("dataHist_B"));

  TCanvas* c_C =  new TCanvas("canvas_C","My plots",1108,4,550,520);
  c_C->cd();
  TPad *pad_C_1 = new TPad("pad_C_1", "pad_C_1", 0, 0.16, 0.98, 1.0);
  pad_C_1->SetTicks(1,1);
  pad_C_1->Draw(); pad_C_1->cd();
  RooPlot* myPlot_C = ws->var("ctau3DRes")->frame(nCtauResBins); // bins
  myPlot_C->SetTitle("");
  //ws->data("dsAB")->plotOn(myPlot_C,Name("dataHist_C"));

  TCanvas* c_D =  new TCanvas("canvas_D","My plots",4,565,550,520);
  c_D->cd();
  TPad *pad_D_1 = new TPad("pad_D_1", "pad_D_1", 0, 0.16, 0.98, 1.0);
  pad_D_1->SetTicks(1,1);
  pad_D_1->Draw(); pad_D_1->cd();
  RooPlot* myPlot_D = wsmc->var("ctau3Dtrue")->frame(Bins(nCtauTrueBins), Range(-1, 6)); // bins
  ws->data("reducedDS_MC")->plotOn(myPlot_D,Name("mcHist_D"));
  myPlot_D->SetTitle("");

  TCanvas* c_E =  new TCanvas("canvas_E","My plots",554,565,550,520);
  c_E->cd();
  TPad *pad_E_1 = new TPad("pad_E_1", "pad_E_1", 0, 0.16, 0.98, 1.0);
  pad_E_1->SetTicks(1,1);
  pad_E_1->Draw(); pad_E_1->cd();
  RooPlot* myPlot_E = ws->var("ctau3D")->frame(Bins(nCtauBins), Range(ctauLow, ctauHigh)); // bins
  myPlot_E->SetTitle("");
  //ws->data("dsAB")->plotOn(myPlot_C,Name("dataHist_C"));

  TCanvas* c_F =  new TCanvas("canvas_F","My plots",1108,565,550,520);
  c_F->cd();
  TPad *pad_F_1 = new TPad("pad_F_1", "pad_F_1", 0, 0.16, 0.98, 1.0);
  pad_F_1->SetTicks(1,1);
  pad_F_1->Draw(); pad_F_1->cd();
  RooPlot* myPlot_F = ws->var("ctau3D")->frame(nCtauBins); // bins
  myPlot_F->SetTitle("");
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

  pad_A_1->cd();
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
  myPlot2_A->GetYaxis()->SetRangeUser(ws->var("nSig")->getVal()/1000, ws->var("nSig")->getVal());
  //myPlot2_A->SetMinimum(2*10);
  myPlot2_A->GetXaxis()->SetLabelSize(0);
  myPlot2_A->GetXaxis()->SetTitleSize(0);
  myPlot2_A->GetXaxis()->CenterTitle();
  myPlot2_A->GetXaxis()->SetRangeUser(massLow,massHigh);
  myPlot2_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  myPlot2_A->Draw();

  float text_x = 0.15;
  float text_y = 0.816;
  float y_diff = 0.05;
  float text_size = 13;
  int text_color = 1;
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);
  drawText(Form("n_{J/#psi} = %.f #pm %.f",ws->var("nSig")->getVal(),ws->var("nSig")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
  drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("nBkg")->getVal(),ws->var("nBkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);

  TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0, 0.006, 0.98, 0.227);
  c_A->cd();
  pad_A_2->Draw();
  pad_A_2->cd();
  pad_A_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_A_2->SetBottomMargin(0.67);
  pad_A_2->SetBottomMargin(0.4);
  pad_A_2->SetFillStyle(4000);
  pad_A_2->SetFrameFillStyle(4000);
  pad_A_2->SetTicks(1,1);

  RooPlot* frameTMP = (RooPlot*)myPlot2_A->Clone("TMP");
  RooHist* hpull_A = frameTMP->pullHist("dataOS_FIT_A","modelHist_A", true);
  hpull_A->SetMarkerSize(0.8);
  RooPlot* pullFrame_A = ws->var("mass")->frame(Title("Pull Distribution")) ;
  pullFrame_A->addPlotable(hpull_A,"P") ;
  pullFrame_A->SetTitle("");
  pullFrame_A->SetTitleSize(0);
  pullFrame_A->GetYaxis()->SetTitleOffset(0.3) ;
  pullFrame_A->GetYaxis()->SetTitle("Pull") ;
  pullFrame_A->GetYaxis()->SetTitleSize(0.15) ;
  pullFrame_A->GetYaxis()->SetLabelSize(0.15) ;
  pullFrame_A->GetYaxis()->SetRangeUser(-3.8,3.8);
  pullFrame_A->GetYaxis()->CenterTitle();

  pullFrame_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame_A->GetXaxis()->SetTitleOffset(1.05) ;
  pullFrame_A->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame_A->GetXaxis()->SetLabelSize(0.15) ;
  pullFrame_A->GetXaxis()->SetTitleSize(0.15) ;
  pullFrame_A->GetXaxis()->CenterTitle();

  pullFrame_A->GetYaxis()->SetTickSize(0.04);
  pullFrame_A->GetYaxis()->SetNdivisions(404);
  pullFrame_A->GetXaxis()->SetTickSize(0.03);
  pullFrame_A->Draw() ;

  TLine *l1 = new TLine(massLow,0,massHigh,0);
  l1->SetLineStyle(1);
  l1->Draw("same");

  printChi2(ws, pad_A_2, frameTMP, fitResSim, "mass", "dataOS_FIT_A", "modelHist_A", nMassBin, false);

  //c_A->Draw();
  //pad_A_1->Draw();
  //pad_A_2->Draw();
  //pad_A_1->Update();
  //pad_A_2->Update();

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
  auto clone_mass_pdf = (RooAbsPdf*)cloneSet->find("model_A");
  clone_mass_pdf->setOperMode(RooAbsArg::ADirty, kTRUE);

  RooStats::SPlot sData = RooStats::SPlot("sData", "An SPlot", *data, clone_mass_pdf, yieldList);
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

  pad_B_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot2_B = (RooPlot*)myPlot_B->Clone();
  ws->data("dataset_SPLOT")->plotOn(myPlot2_B,Name("dataHist_Tot"), MarkerSize(.7), Binning(nCtauErrBins));
  ws->pdf("totPdf")->plotOn(myPlot2_B,Name("totPdf"), LineColor(kGreen+1), Range("ctauErrRange"), LineWidth(2));
  ws->data("dataw_Sig")->plotOn(myPlot2_B,Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2), Binning(nCtauErrBins));
  ws->pdf("sigPdf")->plotOn(myPlot2_B,Name("sigPdf"),LineColor(kRed+2), LineWidth(2), Range("ctauErrRange"));
  ws->data("dataw_Bkg")->plotOn(myPlot2_B,Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue+2), MarkerColor(kBlue+2), Binning(nCtauErrBins));
  ws->pdf("bkgPdf")->plotOn(myPlot2_B,Name("bkgPdf"), LineColor(kBlue+2), LineWidth(2), Range("ctauErrRange"));
  myPlot2_B->GetYaxis()->SetRangeUser(10e-4, 10e8);
  myPlot2_B->GetXaxis()->CenterTitle();
  myPlot2_B->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
  myPlot2_B->SetFillStyle(4000);
  myPlot2_B->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_B->GetXaxis()->SetLabelSize(0);
  myPlot2_B->GetXaxis()->SetTitleSize(0);
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
  TPad *pad_B_2 = new TPad("pad_B_2", "pad_B_2", 0, 0.006, 0.98, 0.227);
  c_B->cd();
  pad_B_2->Draw();
  pad_B_2->cd();
  pad_B_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_B_2->SetBottomMargin(0.67);
  pad_B_2->SetBottomMargin(0.4);
  pad_B_2->SetFillStyle(4000);
  pad_B_2->SetFrameFillStyle(4000);
  pad_B_2->SetTicks(1,1);

  RooPlot* frameTMP_B = (RooPlot*)myPlot2_B->Clone("TMP");
  RooHist* hpull_B = frameTMP_B->pullHist("dataHist_Tot","totPdf");
  hpull_B->SetMarkerSize(0.8);
  RooPlot* pullFrame_B = ws->var("ctau3DErr")->frame(Title("Pull Distribution")) ;
  pullFrame_B->addPlotable(hpull_B,"PX") ;
  pullFrame_B->SetTitle("");
  pullFrame_B->SetTitleSize(0);
  pullFrame_B->GetYaxis()->SetTitleOffset(0.3) ;
  pullFrame_B->GetYaxis()->SetTitle("Pull") ;
  pullFrame_B->GetYaxis()->SetTitleSize(0.15) ;
  pullFrame_B->GetYaxis()->SetLabelSize(0.15) ;
  pullFrame_B->GetYaxis()->SetRangeUser(-3.8,3.8);
  pullFrame_B->GetYaxis()->CenterTitle();

  pullFrame_B->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
  pullFrame_B->GetXaxis()->SetTitleOffset(1.05) ;
  pullFrame_B->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame_B->GetXaxis()->SetLabelSize(0.15) ;
  pullFrame_B->GetXaxis()->SetTitleSize(0.15) ;
  pullFrame_B->GetXaxis()->CenterTitle();

  pullFrame_B->GetYaxis()->SetTickSize(0.04);
  pullFrame_B->GetYaxis()->SetNdivisions(404);
  pullFrame_B->GetXaxis()->SetTickSize(0.03);
  pullFrame_B->Draw() ;

  TLine *lB = new TLine(ctauErrLow,0, ctauErrHigh,0);
  lB->SetLineStyle(1);
  lB->Draw("same");

  cout << endl << "************** Finished SPLOT *****************" << endl << endl;
  //***********************************************************************
  //**************************** CTAU RES FIT *****************************
  //***********************************************************************
  cout << endl << "*************** Start Res Fit *****************" << endl << endl;

  RooDataSet *ctauResCutDS =(RooDataSet*)dataw_Sig->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")),*(ws->var("ctau3DErr"))),"ctau3DRes<0");
  ctauResCutDS->SetName("ctauResCutDS");
  ws->import(*ctauResCutDS);
  double entries = ws->data("ctauResCutDS")->sumEntries();
  //double entries = ws->var("nSig")->getVal();
  cout<<"[Info] #J/psi: "<<entries<<endl;
  ws->factory(Form("N_Jpsi[%.12f,%.12f,%.12f]", entries, entries, entries*2.0));
// create the variables for this model
  //ws->var("ctau3DRes")->setRange("ctauResRange", ctauResLow, ctauResHigh);
  int nGauss = 3;
  ws->factory("ctauRes_mean[0.0]");
  ws->factory("ctau1_CtauRes[0., -0.1, 0.1]");  ws->factory("s1_CtauRes[.5, 0., 3.]");
  ws->factory("ctau2_CtauRes[0., -0.1, 0.1]");  ws->factory("s2_CtauRes[1.26, 0., 2.]");
  ws->factory("ctau3_CtauRes[0., -0.1, 0.1]");  ws->factory("s3_CtauRes[3.37, 0., 7.]");
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
  
  RooAbsPdf *ctauResModel = ctauResModel = new RooAddPdf("GaussModel_Tot", "GaussModel_Tot", *ws->pdf("GaussModelCOND_ctauRes"), *ws->var("N_Jpsi"));
  ws->import(*ctauResModel);

  pad_C_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot2_C = (RooPlot*)myPlot_C->Clone();
  bool isWeighted = ws->data("ctauResCutDS")->isWeighted();
  RooFitResult* fitCtauRes = ws->pdf("GaussModel_Tot")->fitTo(*ctauResCutDS, Save(), Range(ctauResLow,0), SumW2Error(isWeighted), Extended(kTRUE), NumCPU(12), PrintLevel(-1));
  ws->data("ctauResCutDS")->plotOn(myPlot2_C, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0), 
      MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2), Precision(1e-5));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_ctauRes"), Precision(1e-5), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), LineColor(kBlack));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm1"), Precision(1e-5), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel1_ctauRes")), LineColor(kGreen+2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm2"), Precision(1e-5), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel2_ctauRes")), LineColor(kRed+2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm3"), Precision(1e-5), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel3_ctauRes")), LineColor(kBlue+2));
  if(nGauss==4){ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm4"), Precision(1e-5), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel4_ctauRes")), LineColor(kMagenta+2));}
  ws->data("ctauResCutDS")->plotOn(myPlot2_C, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0), 
      MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2), Precision(1e-5));
  myPlot2_C->GetYaxis()->SetRangeUser(10e-4, 10e6);
  myPlot2_C->GetXaxis()->CenterTitle();
  myPlot2_C->GetXaxis()->SetTitle("#frac{l_{J/#psi}}{#sigma_{J/#psi}}");
  myPlot2_C->SetFillStyle(4000);
  myPlot2_C->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_C->GetXaxis()->SetLabelSize(0);
  myPlot2_C->GetXaxis()->SetTitleSize(0);
  myPlot2_C->Draw();

  TLegend* leg_C = new TLegend(text_x+0.29,text_y+0.03,text_x+0.39,text_y-0.17); leg_C->SetTextSize(text_size);
  leg_C->SetTextFont(43);
  leg_C->SetBorderSize(0);
  leg_C->AddEntry(myPlot2_C->findObject("dataHist_ctauRes"),"Data_Sig","pe");
  leg_C->AddEntry(myPlot2_C->findObject("modelHist_ctauRes"),"Total PDF","l");
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
  drawText(Form("s1_{Res} = %.4f #pm %.4f", ws->var("s1_CtauRes")->getVal(), ws->var("s1_CtauRes")->getError() ),text_x+0.5,text_y-y_diff*1,text_color,text_size);
  drawText(Form("s2_{Res} = %.4f #pm %.4f", ws->var("s2_CtauRes")->getVal(), ws->var("s2_CtauRes")->getError() ),text_x+0.5,text_y-y_diff*2,text_color,text_size);
  drawText(Form("s3_{Res} = %.4f #pm %.4f", ws->var("s3_CtauRes")->getVal(), ws->var("s3_CtauRes")->getError() ),text_x+0.5,text_y-y_diff*3,text_color,text_size);
  drawText(Form("f2_{Res} = %.4f #pm %.4f", ws->var("f2_CtauRes")->getVal(), ws->var("f2_CtauRes")->getError() ),text_x+0.5,text_y-y_diff*4,text_color,text_size);
  drawText(Form("f_{Res} = %.4f #pm %.4f", ws->var("f_CtauRes")->getVal(), ws->var("f_CtauRes")->getError() ),text_x+0.5,text_y-y_diff*5,text_color,text_size);
  //drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c", ws->ws->var("s2_CtauRes")->getVal()/ws->var("s1_CtauRes")->getVal(), ptHigh ),text_x+0.5,text_y,text_color,text_size);

  TPad *pad_C_2 = new TPad("pad_C_2", "pad_C_2", 0, 0.006, 0.98, 0.227);
  c_C->cd();
  pad_C_2->Draw();
  pad_C_2->cd();
  pad_C_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_C_2->SetBottomMargin(0.67);
  pad_C_2->SetBottomMargin(0.4);
  pad_C_2->SetFillStyle(4000);
  pad_C_2->SetFrameFillStyle(4000);
  pad_C_2->SetTicks(1,1);

  RooPlot* frameTMP_C = (RooPlot*)myPlot2_C->Clone("TMP");
  RooHist* hpull_C = frameTMP_C->pullHist("dataHist_ctauRes","modelHist_ctauRes", true);
  hpull_C->SetMarkerSize(0.8);
  RooPlot* pullFrame_C = ws->var("ctau3DRes")->frame(Title("Pull Distribution"), Bins(nCtauResBins), Range(ctauResLow, ctauResHigh)) ;
  pullFrame_C->addPlotable(hpull_C,"PX") ;
  pullFrame_C->SetTitle("");
  pullFrame_C->SetTitleSize(0);
  pullFrame_C->GetYaxis()->SetTitleOffset(0.3) ;
  pullFrame_C->GetYaxis()->SetTitle("Pull") ;
  pullFrame_C->GetYaxis()->SetTitleSize(0.15) ;
  pullFrame_C->GetYaxis()->SetLabelSize(0.15) ;
  pullFrame_C->GetYaxis()->SetRangeUser(-3.8,3.8);
  pullFrame_C->GetYaxis()->CenterTitle();

  pullFrame_C->GetXaxis()->SetTitle("#frac{l_{J/#psi}}{#sigma_{J/#psi}}");
  pullFrame_C->GetXaxis()->SetTitleOffset(1.05) ;
  pullFrame_C->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame_C->GetXaxis()->SetLabelSize(0.15) ;
  pullFrame_C->GetXaxis()->SetTitleSize(0.15) ;
  pullFrame_C->GetXaxis()->CenterTitle();

  pullFrame_C->GetYaxis()->SetTickSize(0.04);
  pullFrame_C->GetYaxis()->SetNdivisions(404);
  pullFrame_C->GetXaxis()->SetTickSize(0.03);
  pullFrame_C->Draw() ;

  TLine *lC = new TLine(ctauResLow,0, ctauResHigh,0);
  lC->SetLineStyle(1);
  lC->Draw("same");

  printChi2(ws, pad_C_2, frameTMP_C, fitCtauRes, "ctau3DRes", "dataHist_ctauRes", "modelHist_ctauRes", nCtauResBins, false);
  pad_C_2->Update();
  cout << endl << "************* Finished Sig Res Fit ****************" << endl << endl;
  //***********************************************************************
  //**************************** CTAU TRUE FIT ****************************
  //***********************************************************************
  cout << endl << "************ Start MC Ctau True Fit ***************" << endl << endl;
//MC NP ctau true
  double entries_True = ws->data("reducedDS_MC")->numEntries();
  ws->factory(Form("N_Jpsi_MC[%.12f,%.12f,%.12f]", entries_True, entries_True*0.9, entries));
  ws->factory("lambdaDSS[0.2, 1e-6, 1.]");
  // create the PDF
  ws->factory(Form("TruthModel::%s(%s)", "TruthModel_ctauTrue", "ctau3Dtrue"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "ctauTruePdf", "ctau3Dtrue",  
        "lambdaDSS", 
        "TruthModel_ctauTrue"));

  RooAbsPdf *ctauTrueModel = ctauTrueModel = new RooAddPdf("TrueModel_Tot", "TrueModel_Tot", *ws->pdf("ctauTruePdf"), *ws->var("N_Jpsi_MC"));
  ws->import(*ctauTrueModel);
  
  //c_D->cd();
  //c_D->SetLogy();
  pad_D_1->cd();
  gPad->SetLogy();
  RooFitResult* fitCtauTrue = ws->pdf("TrueModel_Tot")->fitTo(*reducedDS_MC, Save(), Range("ctauTrueRange"), Extended(kTRUE), NumCPU(12), PrintLevel(3));
  RooPlot* myPlot2_D = (RooPlot*)myPlot_D->Clone();
  ws->data("reducedDS_MC")->plotOn(myPlot2_D,Name("MCHist_Tot"), DataError(RooAbsData::SumW2), XErrorSize(0), Precision(1e-5), 
      MarkerSize(.7), Binning(nCtauTrueBins));
  ws->pdf("TrueModel_Tot")->plotOn(myPlot2_D,Name("MCpdf_Tot"), Precision(1e-5), LineColor(kRed+2));
  myPlot2_D->GetYaxis()->SetRangeUser(10e-2, ws->data("reducedDS_MC")->sumEntries()*10);
  myPlot2_D->GetXaxis()->SetRangeUser(-1, 7);
  myPlot2_D->GetXaxis()->CenterTitle();
  myPlot2_D->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} MC True (mm)");
  myPlot2_D->SetFillStyle(4000);
  myPlot2_D->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_D->GetXaxis()->SetLabelSize(0);
  myPlot2_D->GetXaxis()->SetTitleSize(0);
  myPlot2_D->Draw();
  
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),text_x+0.5,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x+0.5,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x+0.5,text_y-y_diff,text_color,text_size);
  drawText(Form("#lambdaDSS = %.4f #pm %.4f", ws->var("lambdaDSS")->getVal(), ws->var("lambdaDSS")->getError() ),text_x+0.5,text_y-y_diff*2,text_color,text_size);

  //TPad *pad_D_2 = new TPad("pad_D_2", "pad_D_2", 0, 0.006, 0.98, 0.227);
  //RooPlot* frameTMP_D = (RooPlot*)myPlot2_D->Clone("TMP");
  //RooHist* hpull_D;
  //pullDist(ws, pad_D_2, c_D, frameTMP_D, hpull_D, "MCHist_Tot", "MCpdf_Tot", "ctau3Dtrue", nCtauTrueBins, -1, 6, "#font[12]{l}_{J/#psi} MC True (mm)");
  //printChi2_test(ws, pad_D_2, frameTMP_D, hpull_D, fitCtauTrue, "ctau3Dtrue", "MCHist_Tot", "MCpdf_Tot", nCtauTrueBins);

  TPad *pad_D_2 = new TPad("pad_D_2", "pad_D_2", 0, 0.006, 0.98, 0.227);
  c_D->cd();
  pad_D_2->Draw();
  pad_D_2->cd();
  pad_D_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_D_2->SetBottomMargin(0.67);
  pad_D_2->SetBottomMargin(0.4);
  pad_D_2->SetFillStyle(4000);
  pad_D_2->SetFrameFillStyle(4000);
  pad_D_2->SetTicks(1,1);

  RooPlot* frameTMP_D = (RooPlot*)myPlot2_D->Clone("TMP");
  RooHist* hpull_D = frameTMP_D->pullHist(0,0,true);
  hpull_D->SetMarkerSize(0.8);
  RooPlot* pullFrame_D = ws->var("ctau3Dtrue")->frame(Title("Pull Distribution"), Bins(nCtauTrueBins), Range(-1, 6.5)) ;
  pullFrame_D->addPlotable(hpull_D,"PX") ;
  pullFrame_D->SetTitle("");
  pullFrame_D->SetTitleSize(0);
  pullFrame_D->GetYaxis()->SetTitleOffset(0.3) ;
  pullFrame_D->GetYaxis()->SetTitle("Pull") ;
  pullFrame_D->GetYaxis()->SetTitleSize(0.15) ;
  pullFrame_D->GetYaxis()->SetLabelSize(0.15) ;
  pullFrame_D->GetYaxis()->SetRangeUser(-3.8,3.8);
  pullFrame_D->GetYaxis()->CenterTitle();

  pullFrame_D->GetXaxis()->SetTitle("#frac{l_{J/#psi}}{#sigma_{J/#psi}}");
  pullFrame_D->GetXaxis()->SetTitleOffset(1.05) ;
  pullFrame_D->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame_D->GetXaxis()->SetLabelSize(0.15) ;
  pullFrame_D->GetXaxis()->SetTitleSize(0.15) ;
  pullFrame_D->GetXaxis()->CenterTitle();

  pullFrame_D->GetYaxis()->SetTickSize(0.04);
  pullFrame_D->GetYaxis()->SetNdivisions(404);
  pullFrame_D->GetXaxis()->SetTickSize(0.03);
  pullFrame_D->Draw() ;

  TLine *lD = new TLine(-1, 0, 6.5, 0);
  lD->SetLineStyle(1);
  lD->Draw("same");

  printChi2(ws, pad_D_2, frameTMP_D, fitCtauTrue, "ctau3Dtrue", "MCHist_Tot", "MCpdf_Tot", nCtauTrueBins, false);
  pad_D_2->Update();
  //***********************************************************************
  //**************************** Bkg CTAU FIT *****************************
  //***********************************************************************
  cout << endl << "************** Start Ctau Fit *****************" << endl << endl;
//make parameter 3 exp
  ws->factory("zeroMean[0.]");
  ws->factory("fDFSS[0.5, 1e-6, 1.]");
  ws->factory("fDLIV[0.5, 1e-6, 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg[0.3, 0.3, 1.]");
  ws->factory("lambdaDSS_Bkg[0.5, 1e-6, 1.]");
//parameters fixed by Resolution model
  ws->var("ctauRes_mean")->setConstant(kTRUE);
  ws->var("s1_CtauRes")->setConstant(kTRUE);
  ws->var("s2_CtauRes")->setConstant(kTRUE);
  ws->var("s3_CtauRes")->setConstant(kTRUE);
  ws->var("f2_CtauRes")->setConstant(kTRUE);
  ws->var("f_CtauRes")->setConstant(kTRUE); 
  ws->var("zeroMean")->setConstant(kTRUE);
//make res model
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
  ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "pdfCTAURES", "ctauRes32", "ctauRes1", "f_CtauRes"));
//make 3 exp
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS", "ctau3D",
                 "lambdaDSS_Bkg", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF", "ctau3D",
                 "lambdaDF_Bkg", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", "pdfCTAUDDS", "ctau3D",
                  "lambdaDDS_Bkg", "pdfCTAURES"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU1",
                  "fDFSS", "pdfCTAUDSS", "pdfCTAUDF"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfBkgCTAUCONDNoPR",//NP
                  "fDLIV", "pdfCTAU1", "pdfCTAUDDS"));
  ws->factory(Form("SUM::%s(%s)", "pdfBkgCTAUCONDPR",//PR
        "pdfCTAURES"//resolution model
        ));
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//make res model
  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauResSig1", "ctau3D", 
                   "ctauRes_mean", //"ctau1_CtauRes",
                   "s1_CtauRes",
                   "zeroMean",
                   "ctau3DErr"
                   ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauResSig2", "ctau3D", 
                   "ctauRes_mean", //"ctau2_CtauRes",
                   "s2_CtauRes",
                   "zeroMean",
                   "ctau3DErr"
                   ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauResSig3", "ctau3D", 
                   "ctauRes_mean", //"ctau3_CtauRes",
                   "s3_CtauRes",
                   "zeroMean",
                   "ctau3DErr"
                   ));
  ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "ctauResSig32", "ctauResSig3", "ctauResSig2", "f2_CtauRes"));
  ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "pdfSigCTAURES", "ctauResSig32", "ctauResSig1", "f_CtauRes"));
//make parameter 3 exp
  ws->var("lambdaDSS")->setConstant(kTRUE);
  ws->factory("fDFSS_Sig[0.5, 1e-6, 1.]");
  ws->factory("fDLIV_Sig[0.5, 1e-6, 1.]");
  ws->factory("lambdaDDS_Sig[0.03, 0.01, 1.]");
  ws->factory("lambdaDF_Sig[0.04, 0.01, 1.]");
  //ws->factory("lambdaDSS_Sig[0.5, 1e-6, 1.]");
//make 3 exp
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfSigCTAUDSS", "ctau3D",
                 "lambdaDSS", "pdfSigCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfSigCTAUDF", "ctau3D",
                 "lambdaDF_Sig", "pdfSigCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", "pdfSigCTAUDDS", "ctau3D",
                  "lambdaDDS_Sig", "pdfSigCTAURES"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfSigCTAU1",
                  "fDFSS_Sig", "pdfSigCTAUDSS", "pdfSigCTAUDF"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfSigCTAUCONDNoPR",//NP
                  "fDLIV_Sig", "pdfSigCTAU1", "pdfSigCTAUDDS"));
  ws->factory(Form("SUM::%s(%s)", "pdfSigCTAUCONDPR",//PR
                  "pdfSigCTAURES"//resolution model
                  ));
//multiply to mass pdf
//  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgNoPR", "pdfBkgCTAUCONDNoPR", "bkg" ));//bkg is background mass pdf.
//  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgPR", "pdfBkgCTAUCONDPR", "bkg" ));
////??
//  RooProdPdf pdfNoPR("pdfCTAUMASS_BkgNoPR", "", 
//      *ws->pdf("bkgPdf"),//ctau Err template
//      Conditional( *ws->pdf("pdfBkgCTAUCONDNoPR"), RooArgList(*ws->var("ctau3D")) )
//      );
//  RooProdPdf pdfPR("pdfCTAUMASS_BkgPR", "", 
//      *ws->pdf("bkgPdf"), 
//      Conditional( *ws->pdf("pdfBkgCTAUCONDPR"), RooArgList(*ws->var("ctau3D")) )
//      );
//  ws->import(pdfNoPR); ws->import(pdfPR); //import
//b_Bkg fraction  
  ws->factory("b_Bkg[0.6, 0., 1.]");//NP fraction for bkg
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfBkgCTAUCOND",
                  "b_Bkg",
                  "pdfBkgCTAUCONDNoPR",
                  "pdfBkgCTAUCONDPR"
                  ));
  ws->factory(Form("RooExtendPdf::%s(%s,%s)", "pdfTot_Bkg", "pdfBkgCTAUCOND", "nBkg"));//nBkg is number of bkg from dataw_Bkg

  ws->factory("b_Frac[0.6, 0., 1.]");//NP fraction for Sig
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfSigCTAUCOND",
                  "b_Frac",
                  "pdfSigCTAUCONDNoPR",
                  "pdfSigCTAUCONDPR"
                  ));
  ws->factory(Form("RooExtendPdf::%s(%s,%s)", "pdfTot_Sig", "pdfSigCTAUCOND", "nSig"));//nSig is number of bkg from dataw_Sig
//multiply to mass pdf
//  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgNoPR", "pdfBkgCTAUCONDNoPR", "bkg" ));//bkg is background mass pdf.
//  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgPR", "pdfBkgCTAUCONDPR", "bkg" ));
////??
//  RooProdPdf pdfNoPR("pdfCTAUMASS_BkgNoPR", "", 
//      *ws->pdf("bkgPdf"),//ctau Err template
//      Conditional( *ws->pdf("pdfBkgCTAUCONDNoPR"), RooArgList(*ws->var("ctau3D")) )
//      );
//  RooProdPdf pdfPR("pdfCTAUMASS_BkgPR", "", 
//      *ws->pdf("bkgPdf"), 
//      Conditional( *ws->pdf("pdfBkgCTAUCONDPR"), RooArgList(*ws->var("ctau3D")) )
//      );
//  ws->import(pdfNoPR); ws->import(pdfPR); //import
//b_Bkg fraction  
  //c_E->cd();
  //c_E->SetLogy();
  pad_E_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot2_E = (RooPlot*)myPlot_E->Clone();
  RooFitResult* fitCtauBkg = ws->pdf("pdfTot_Bkg")->fitTo(*dataw_Bkg, Save(), Range("ctauRange"), Extended(kTRUE), NumCPU(12), PrintLevel(-1));
  //ws->data("dataw_Sig")->plotOn(myPlot2_E,Name("dataHist_Tot"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kRed+2), MarkerColor(kRed+2));
  ws->data("dataw_Bkg")->plotOn(myPlot2_E,Name("dataHist_ctauBkg"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kBlue+2), MarkerColor(kBlue+2));
  //ws->pdf("BkgModel_Tot")->plotOn(myPlot2_E,Name("ctauBkg_Tot"), LineColor(kBlack));
  ws->pdf("pdfTot_Bkg")->plotOn(myPlot2_E,Name("ctauBkg_Tot"), 
      Normalization(ws->data("dataw_Bkg")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE),
      FillStyle(1001), FillColor(kAzure-9), LineColor(kBlack), LineWidth(3), DrawOption("LCF"), Precision(1e-4));
      //FillStyle(1001), FillColor(kAzure-9), LineColor(kBlack), VLines(), DrawOption("LCF"), Precision(1e-4));
  ////ws->pdf("BkgModel_Tot")->plotOn(myPlot2_E,Name("BKG"), 
  ////    Normalization(ws->data("dataw_Bkg")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
  ////    ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE), Components("pdfCTAU1"),
  ////    LineWidth(2), LineColor(kRed+2));
  ws->pdf("pdfTot_Bkg")->plotOn(myPlot2_E,Name("BKG"), 
      Normalization(ws->data("dataw_Bkg")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE), Components("pdfBkgCTAUCONDPR"),
      LineWidth(2), LineColor(kRed+2));
  //ws->pdf("pdfTot_Bkg")->plotOn(myPlot2_E,Name("BKG"), 
  //    Normalization(ws->data("dataw_Bkg")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
  //    ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE), Components("pdfCTAUMASS_BkgNoPR"),
  //    LineWidth(2), LineColor(kRed+2));
  ws->data("dataw_Bkg")->plotOn(myPlot2_E,Name("dataHist_ctauBkg"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kBlue+2), MarkerColor(kBlue+2));
  myPlot2_E->GetYaxis()->SetRangeUser(10e-2, 10e6);
  myPlot2_E->GetXaxis()->SetRangeUser(-4, 6);
  myPlot2_E->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  myPlot2_E->SetFillStyle(4000);
  myPlot2_E->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_E->GetXaxis()->SetLabelSize(0);
  myPlot2_E->GetXaxis()->SetTitleSize(0);
  myPlot2_E->Draw();
  TLegend* leg_E = new TLegend(text_x+0.25,text_y+0.03,text_x+0.38,text_y-0.17); leg_E->SetTextSize(text_size);
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
  
  drawText(Form("N_{Bkg} = %.f #pm %.f", ws->var("nBkg")->getVal(), ws->var("nBkg")->getError() ),text_x+0.5,text_y,text_color,text_size);
  drawText(Form("b_{Bkg} = %.4f #pm %.4f", ws->var("b_Bkg")->getVal(), ws->var("b_Bkg")->getError() ),text_x+0.5,text_y-y_diff*1,text_color,text_size);
  drawText(Form("fDFSS = %.4f #pm %.4f", ws->var("fDFSS")->getVal(), ws->var("fDFSS")->getError() ),text_x+0.5,text_y-y_diff*2,text_color,text_size);
  drawText(Form("fDLIV = %.4f #pm %.4f", ws->var("fDLIV")->getVal(), ws->var("fDLIV")->getError() ),text_x+0.5,text_y-y_diff*3,text_color,text_size);
  drawText(Form("#lambdaDDS_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDDS_Bkg")->getVal(), ws->var("lambdaDDS_Bkg")->getError() ),text_x+0.5,text_y-y_diff*4,text_color,text_size);
  drawText(Form("#lambdaDF_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDF_Bkg")->getVal(), ws->var("lambdaDF_Bkg")->getError() ),text_x+0.5,text_y-y_diff*5,text_color,text_size);
  drawText(Form("#lambdaDSS_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDSS_Bkg")->getVal(), ws->var("lambdaDSS_Bkg")->getError() ),text_x+0.5,text_y-y_diff*6,text_color,text_size);

  TPad *pad_E_2 = new TPad("pad_E_2", "pad_E_2", 0, 0.006, 0.98, 0.227);
  RooPlot* frameTMP_E = (RooPlot*)myPlot2_E->Clone("TMP");
  RooHist* hpull_E;
  pullDist(ws, pad_E_2, c_E, frameTMP_E, hpull_E, "dataHist_ctauBkg", "ctauBkg_Tot", "ctau3D", nCtauBins, ctauLow, ctauHigh, "#font[12]{l}_{J/#psi} (mm)");
  printChi2_test(ws, pad_E_2, frameTMP_E, hpull_E, fitCtauBkg, "ctau3D", "dataHist_ctauBkg", "ctauBkg_Tot", nCtauBins);
  pad_E_2->Update();
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
  //ws->factory("b_Bkg[0.5, 1e-6, 1.0]");
  //ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU_Bkg",
  //                "b_Bkg",
  //                "pdfCTAU_BkgNoPR",
  //                "pdfCTAU_BkgPR"
  //                ));
  //ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfBkgCTAUCOND",
  //                "b_Bkg",
  //                "pdfBkgCTAUCONDNoPR",
  //                "pdfBkgCTAUCONDPR"
  //                ));
  //ws.factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgPR",
  //                "pdfBkgCTAUCONDPR",
  //                "bkg"
  //                ));
  //ws.factory(Form("RooExtendPdf::%s(%s,%s)", "pdfCTAUMASS_Tot_BkgPR",
  //                "pdfCTAUMASS_BkgPR",
  //                "nBkg"
  //                ));
  //ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_JpsiPR",
  //               Form("%s_JpsiPR", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"),
  //               "pdfMASS_Jpsi_%s",
  //               ));

  //c_F->cd();
  //c_F->SetLogy();
  pad_F_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot2_F = (RooPlot*)myPlot_F->Clone();
  RooFitResult* fitCtauSig = ws->pdf("pdfTot_Sig")->fitTo(*dataw_Sig, Save(), Range("ctauRange"), Extended(kTRUE), NumCPU(12), PrintLevel(-1));
  ws->data("dataw_Sig")->plotOn(myPlot2_F,Name("dataHist_ctauSig"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kRed+2), MarkerColor(kRed+2));
  //ws->pdf("SigModel_Tot")->plotOn(myPlot_F,Name("ctauSig_Tot"), LineColor(kBlack));
  ws->pdf("pdfTot_Sig")->plotOn(myPlot2_F, Name("ctauSig_Tot"), 
      Normalization(ws->data("dataw_Sig")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Sig"), kTRUE),
      FillStyle(1001), FillColor(kRed-10), LineColor(kBlack), LineWidth(3), DrawOption("LCF"), Precision(1e-4));
      //LineColor(kBlack), LineWidth(3), DrawOption("LCF"), Precision(1e-4));
      //FillStyle(1001), FillColor(kAzure-9), LineColor(kBlack), VLines(), DrawOption("LCF"), Precision(1e-4));
  ////ws->pdf("SigModel_Tot")->plotOn(myPlot_F,Name("BKG"), 
  ////    Normalization(ws->data("dataw_Sig")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
  ////    ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Sig"), kTRUE), Components("pdfCTAU1"),
  ////    LineWidth(2), LineColor(kRed+2));
  ws->pdf("pdfTot_Sig")->plotOn(myPlot2_F,Name("PR"), 
      Normalization(ws->data("dataw_Sig")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Sig"), kTRUE), Components("pdfSigCTAUCONDPR"),
      LineWidth(2), LineColor(kMagenta+3));
  ws->pdf("pdfTot_Sig")->plotOn(myPlot2_F,Name("NP"), 
      Normalization(ws->data("dataw_Sig")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Sig"), kTRUE), Components("pdfSigCTAUCONDNoPR"),
      LineWidth(2), LineColor(kOrange-5));
  //ws->pdf("pdfTot_Sig")->plotOn(myPlot_F,Name("BKG"), 
  //    Normalization(ws->data("dataw_Sig")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
  //    ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Sig"), kTRUE), Components("pdfCTAUMASS_SigNoPR"),
  //    LineWidth(2), LineColor(kRed+2));
  ws->data("dataw_Sig")->plotOn(myPlot2_F,Name("dataHist_ctauSig"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kRed+2), MarkerColor(kRed+2));
  myPlot2_F->GetYaxis()->SetRangeUser(10e-2, 10e6);
  myPlot2_F->GetXaxis()->SetRangeUser(-4, 6);
  myPlot2_F->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  myPlot2_F->SetFillStyle(4000);
  myPlot2_F->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_F->GetXaxis()->SetLabelSize(0);
  myPlot2_F->GetXaxis()->SetTitleSize(0);
  myPlot2_F->Draw();

  TLegend* leg_F = new TLegend(text_x+0.25,text_y+0.03,text_x+0.37,text_y-0.17); leg_F->SetTextSize(text_size);
  leg_F->SetTextFont(43);
  leg_F->SetBorderSize(0);
  leg_F->AddEntry(myPlot2_F->findObject("dataHist_ctauSig"),"Data_Sig","pe");
  leg_F->AddEntry(myPlot2_F->findObject("ctauSig_Tot"),"Total fit","fl");
  leg_F->AddEntry(myPlot2_F->findObject("PR"),"Prompt J/#psi","l");
  leg_F->AddEntry(myPlot2_F->findObject("NP"),"Non-Prompt J/#psi","l");
  //leg_F->AddEntry(myPlot_F->findObject("test"),"?? PDF","l");
  leg_F->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);

  drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("nSig")->getVal(), ws->var("nSig")->getError() ),text_x+0.5,text_y,text_color,text_size);
  drawText(Form("b_{J/#psi} = %.4f #pm %.4f", ws->var("b_Frac")->getVal(),ws->var("b_Frac")->getError() ),text_x+0.5,text_y-y_diff*1,text_color,text_size);
  drawText(Form("fDFSS = %.4f #pm %.4f", ws->var("fDFSS_Sig")->getVal(), ws->var("fDFSS_Sig")->getError() ),text_x+0.5,text_y-y_diff*2,text_color,text_size);
  drawText(Form("fDLIV = %.4f #pm %.4f", ws->var("fDLIV_Sig")->getVal(), ws->var("fDLIV_Sig")->getError() ),text_x+0.5,text_y-y_diff*3,text_color,text_size);
  drawText(Form("#lambdaDDS_{J/#psi} = %.4f #pm %.4f", ws->var("lambdaDDS_Sig")->getVal(), ws->var("lambdaDDS_Sig")->getError() ),text_x+0.5,text_y-y_diff*4,text_color,text_size);
  drawText(Form("#lambdaDF_{J/#psi} = %.4f #pm %.4f", ws->var("lambdaDF_Sig")->getVal(), ws->var("lambdaDF_Sig")->getError() ),text_x+0.5,text_y-y_diff*5,text_color,text_size);
  drawText(Form("#lambdaDSS_{J/#psi} = %.4f #pm %.4f", ws->var("lambdaDSS")->getVal(), ws->var("lambdaDSS")->getError() ),text_x+0.5,text_y-y_diff*6,text_color,text_size);

  TPad *pad_F_2 = new TPad("pad_F_2", "pad_F_2", 0, 0.006, 0.98, 0.227);
  RooPlot* frameTMP_F = (RooPlot*)myPlot2_F->Clone("TMP");
  RooHist* hpull_F;
  pullDist(ws, pad_F_2, c_F, frameTMP_F, hpull_F, "dataHist_ctauSig", "ctauSig_Tot", "ctau3D", nCtauBins, ctauLow, ctauHigh, "#font[12]{l}_{J/#psi} (mm)");
  printChi2_test(ws, pad_F_2, frameTMP_F, hpull_F, fitCtauSig, "ctau3D", "dataHist_ctauSig", "ctauSig_Tot", nCtauBins);
  pad_F_2->Update();
  //***********************************************************************
  //*************************** DRAW CTAU FIT *****************************
  //***********************************************************************
//    c_F->cd();
//    c_F->SetLogy();
//    double normBkg = ws->data("dsAB")->sumEntries()/ws->data("dataw_Bkg")->sumEntries();
//    cout<<"[Info]: "<<normBkg<<endl;
//    RooPlot* myPlot2_F = (RooPlot*)myPlot_F->Clone();
//    myPlot2_F->updateNormVars(RooArgSet(*ws->var("mass"), *ws->var("ctau3D"), *ws->var("ctau3DErr")));
//    ws->pdf("pdfTot_Bkg")->plotOn(myPlot2_F,Name("BKG"), 
//        //Normalization(normBkg, RooAbsReal::NumEvent), NormRange("ctauRange"),
//        Normalization(ws->data("dataw_Bkg")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
//        ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE),
//        FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("CF"), Precision(1e-4));
//    ws->data("dsAB")->plotOn(myPlot2_F,Name("dataHist_ctauTot"), DataError(RooAbsData::SumW2), XErrorSize(0));
//    //ws->pdf("BkgModel_Tot")->plotOn(myPlot2_F,Name("resPdf"), Components(*ws->pdf("pdfCTAUDSS")), LineColor(kRed+2));
//    //ws->pdf("BkgModel_Tot")->plotOn(myPlot2_F,Name("resPdf"), Components(*ws->pdf("pdfCTAUDF")), LineColor(kBlue+2));
//    //ws->pdf("BkgModel_Tot")->plotOn(myPlot2_F,Name("resPdf"), Components(*ws->pdf("pdfCTAUDDS")), LineColor(kMagenta+2));
//    myPlot2_F->GetYaxis()->SetRangeUser(10e-2, 10e6);
//    myPlot2_F->GetXaxis()->SetRangeUser(-4, 6);
//    myPlot2_F->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
//    myPlot2_F->Draw();
//    TLegend* leg_F = new TLegend(text_x+0.5,text_y-0.17,text_x+0.65,text_y-0.03); leg_F->SetTextSize(text_size);
//    leg_F->SetTextFont(43);
//    leg_F->SetBorderSize(0);
//    leg_F->AddEntry(myPlot2_F->findObject("dataHist_ctauTot"),"Data","pe");
//    leg_F->AddEntry(myPlot2_F->findObject("BKG"),"Background PDF","fl");
//    //leg_F->AddEntry(myPlot2_F->findObject("test"),"?? PDF","l");
//    leg_F->Draw("same");
//    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
//    if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
//    else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
//    drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);

  c_A->Update();
  c_B->Update();
  c_C->Update();
  c_D->Update();
  c_E->Update();
  c_F->Update();

  TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

  c_A->SaveAs(Form("figs/2Dfit/Mass_%s.pdf",kineLabel.Data()));
  c_B->SaveAs(Form("figs/2Dfit/ctauErr_%s.pdf",kineLabel.Data()));
  c_C->SaveAs(Form("figs/2Dfit/ctauRes_%s.pdf",kineLabel.Data()));
  c_D->SaveAs(Form("figs/2Dfit/ctauTrue_%s.pdf",kineLabel.Data()));
  c_E->SaveAs(Form("figs/2Dfit/Bkg_%s.pdf",kineLabel.Data()));
  c_F->SaveAs(Form("figs/2Dfit/Sig_%s.pdf",kineLabel.Data()));

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
void printChi2(RooWorkspace* myws, TPad* Pad, RooPlot* frame, RooFitResult* fitRes, string varLabel, string dataLabel, string pdfLabel, int nBins, bool useDefaultName=true)
{
  double chi2=0; unsigned int ndof=0;
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.1);
  //unsigned int nFitPar = myws->pdf(pdfLabel.c_str())->getParameters(*myws->data(dataLabel.c_str()))->selectByAttrib("Constant",kFALSE)->getSize();
  unsigned int nFitPar = fitRes->floatParsFinal().getSize();
  RooHist *hpull = frame->pullHist(dataLabel.c_str(), pdfLabel.c_str(), true);
  double* ypulls = hpull->GetY();
  unsigned int nFullBins = 0;
  for (int i = 0; i < nBins; i++) {
    if(ypulls[i]==0) continue;
      chi2 += ypulls[i]*ypulls[i];
      nFullBins++;
  }
  ndof = nFullBins - nFitPar;
  t->DrawLatex(0.7, 0.85, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));

  if (useDefaultName) {
    RooRealVar chi2Var("chi2","chi2",chi2);
    RooRealVar ndofVar("ndof","ndof",ndof);
    myws->import(chi2Var); myws->import(ndofVar);
  } else {
    RooRealVar chi2Var((string("chi2_")+pdfLabel).c_str(),(string("chi2_")+pdfLabel).c_str(),chi2);
    RooRealVar ndofVar((string("ndof_")+pdfLabel).c_str(),(string("ndof_")+pdfLabel).c_str(),ndof);
    myws->import(chi2Var); myws->import(ndofVar);
  }
  delete hpull;
};

void printChi2_test(RooWorkspace* myws, TPad* Pad, RooPlot* frame, RooHist* hpull, RooFitResult* fitRes, string varLabel, string dataLabel, string pdfLabel, int nBins, bool useDefaultName=true)
{
  double chi2=0; unsigned int ndof=0;
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.1);
  //unsigned int nFitPar = myws->pdf(pdfLabel.c_str())->getParameters(*myws->data(dataLabel.c_str()))->selectByAttrib("Constant",kFALSE)->getSize();
  unsigned int nFitPar = fitRes->floatParsFinal().getSize();
  //RooHist *hpull = frame->pullHist(dataLabel.c_str(), pdfLabel.c_str(), true);
  hpull = frame->pullHist(dataLabel.c_str(), pdfLabel.c_str());
  double* ypulls = hpull->GetY();
  unsigned int nFullBins = 0;
  for (int i = 0; i < nBins; i++) {
    if(ypulls[i]==0) continue;
      chi2 += ypulls[i]*ypulls[i];
      nFullBins++;
  }
  ndof = nFullBins - nFitPar;
  t->DrawLatex(0.7, 0.85, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));

  if (useDefaultName) {
    RooRealVar chi2Var("chi2","chi2",chi2);
    RooRealVar ndofVar("ndof","ndof",ndof);
    myws->import(chi2Var); myws->import(ndofVar);
  } else {
    RooRealVar chi2Var((string("chi2_")+pdfLabel).c_str(),(string("chi2_")+pdfLabel).c_str(),chi2);
    RooRealVar ndofVar((string("ndof_")+pdfLabel).c_str(),(string("ndof_")+pdfLabel).c_str(),ndof);
    myws->import(chi2Var); myws->import(ndofVar);
  }
  delete hpull;
};

void pullDist(RooWorkspace* ws, TPad* Pad, TCanvas* canvas, RooPlot* frame, RooHist* hpull, string dataHist, string modelHist, string variable, int nBins, float Low, float High, string titleX){
  //TPad *Pad = new TPad("Pad", "Pad", 0, 0.006, 0.98, 0.227);
  canvas->cd();
  Pad->Draw();
  Pad->cd();
  Pad->SetTopMargin(0); // Upper and lower plot are joined
  Pad->SetBottomMargin(0.67);
  Pad->SetBottomMargin(0.4);
  Pad->SetFillStyle(4000);
  Pad->SetFrameFillStyle(4000);
  Pad->SetTicks(1,1);

  //RooPlot* frame = (RooPlot*)myPlot2_C->Clone("TMP");
  hpull = frame->pullHist(dataHist.c_str(),modelHist.c_str(), true);
  hpull->SetMarkerSize(0.8);
  RooPlot* pullFrame = ws->var(variable.c_str())->frame(Title("Pull Distribution"), Bins(nBins), Range(Low, High)) ;
  pullFrame->addPlotable(hpull,"PX") ;
  pullFrame->SetTitle("");
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.3) ;
  pullFrame->GetYaxis()->SetTitle("Pull") ;
  pullFrame->GetYaxis()->SetTitleSize(0.15) ;
  pullFrame->GetYaxis()->SetLabelSize(0.15) ;
  pullFrame->GetYaxis()->SetRangeUser(-3.8,3.8);
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle(titleX.c_str());
  pullFrame->GetXaxis()->SetTitleOffset(1.05) ;
  pullFrame->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame->GetXaxis()->SetLabelSize(0.15) ;
  pullFrame->GetXaxis()->SetTitleSize(0.15) ;
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw() ;

  TLine *l = new TLine(Low,0, High,0);
  l->SetLineStyle(1);
  l->Draw("same");

  Pad->Update();
};
