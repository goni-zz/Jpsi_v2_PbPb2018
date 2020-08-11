//This code fits the Jpsi data with either the nominal fit or an alternative fit. The difference between the two fits is the signal shape. The nominal fit fits the signals with double CB functions, while the alternative fit fits them with just a gaussian.

#include <iostream>
#include "rootFitHeaders.h"
#include "commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "cutsAndBin.h"
//#include "PsetCollection.h"
#include "CMS_lumi_v2mass.C"
#include "tdrstyle.C"
//#include "StyleSetting.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"

using namespace std;
using namespace RooFit;
void doFit_Sim_Data( 
       float ptLow=6.5, float ptHigh=50, 
       float yLow=0., float yHigh=0.6,//Run 1 has p going in -z direction
       int cLow=0, int cHigh=200,
       float muPtCut=0.0,
       bool whichModel=0,   // Nominal = 0. Alternative = 1.
       int ICset = 1
			) 
{

  gStyle->SetEndErrorSize(0);

  float massLow = 2.6;
  float massHigh = 3.5;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = 36; //(massHigh-massLow)*30;

  float ctauLow = -3;
  float ctauHigh = 5;

  float ctauLowForPlot = ctauLow;    
  float ctauHighForPlot = ctauHigh;

  int   nCtauBin  = (ctauHigh-ctauLow)*10;

  TFile* f1;

  //The order is {sigma_1,x,alpha_1,n_1,f,err_mu,err_sigma,m_lambda}
  double paramsupper[8] = {0.2, 3.0, 3.321, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  TString kineCut;
  TString SigCut;
  TString BkgCut;

  //Select Data Set
    f1 = new TFile("skimmedFiles/OniaRooDataSet_isMC0_JPsi.root");
    kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f",ptLow, ptHigh, yLow, yHigh);
    //SigCut = Form("pt>%.2f && pt<%.2f && abs(y)<%.2f && abs(eta1)<%.2f && abs(eta2)<%.2f && mass>2.9&&mass<3.2", ptLow, ptHigh, yHigh, yHigh, yHigh);
    //BkgCut = Form("pt>%.2f && pt<%.2f && abs(y)<%.2f && abs(eta1)<%.2f && abs(eta2)<%.2f && (mass<2.9||mass>3.2)",ptLow, ptHigh, yHigh, yHigh, yHigh);
    //TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )";//2018 acceptance cut

  //kineCut = accCut + kineCut;

  //import and merge datasets
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  cout << "####################################" << endl;
  RooDataSet *reducedDS_A = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  //RooDataSet *reducedDS_A = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("ctau3D")), *(ws->var("pt")), *(ws->var("y"))), SigCut.Data() );
  reducedDS_A->SetName("reducedDS_A");
  ws->import(*reducedDS_A);
  RooDataSet *reducedDS_B = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3D")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  //RooDataSet *reducedDS_B = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3D")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), BkgCut.Data() );
  //RooDataSet *reducedDS_B = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3D")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  reducedDS_B->SetName("reducedDS_B");
  ws->import(*reducedDS_B);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();
  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  ws->var("ctau3D")->Print();

  RooCategory tp("tp","tp");
  tp.defineType("A");
  tp.defineType("B");

  // Create a dataset that imports contents of all the above datasets mapped by index category tp
  RooDataSet* dsAB = new RooDataSet("dsAB","dsAB",RooArgSet(*(ws->var("mass")), *(ws->var("ctau3D")), *(ws->var("pt")), *(ws->var("y"))),Index(tp),Import("A",*reducedDS_A),Import("B",*reducedDS_B));
  cout << "******** New Combined Dataset ***********" << endl;
  dsAB->Print();
  ws->import(*dsAB);

  TCanvas* c_A =  new TCanvas("canvas_A","My plots",4,45,550,520);
  c_A->cd();
  RooPlot* myPlot_A = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS_A")->plotOn(myPlot_A,Name("dataHist_A"));

  RooRealVar mean("m_{J/#Psi}","mean of the signal gaussian mass PDF",pdgMass.JPsi, pdgMass.JPsi -0.1, pdgMass.JPsi + 0.1 ) ;

  TCanvas* c_B =  new TCanvas("canvas_B","My plots",504,45,550,520);
  c_B->cd();
  RooPlot* myPlot_B = ws->var("ctau3D")->frame(nCtauBin); // bins
  //ws->data("reducedDS_B")->plotOn(myPlot_B,Name("dataHist_B"));
  RooRealVar mean_B("m_{J/#Psi}_B","mean of the signal gaussian mass PDF",pdgMass.JPsi, pdgMass.JPsi -0.1, pdgMass.JPsi + 0.1 ) ;

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

//A:Mass double CB function
  RooRealVar    sigma_1_A("sigma_1_A","width/sigma of the signal gaussian mass PDF",sigma_1_init, paramslower[0], paramsupper[0]);
  RooRealVar   *x_A = new RooRealVar("x_A","sigma ratio ", x_init, paramslower[1], paramsupper[1]);
  RooFormulaVar sigma_2_A("sigma_2_A","@0*@1",RooArgList(sigma_1_A, *x_A) );
  RooRealVar    alpha_1_A("alpha_1_A","tail shift", alpha_1_init, paramslower[2], paramsupper[2]);
  RooFormulaVar alpha_2_A("alpha_2_A","1.0*@0",RooArgList(alpha_1_A) );
  RooRealVar    n_1_A("n_1_A","power order", n_1_init , paramslower[3], paramsupper[3]);
  RooFormulaVar n_2_A("n_2_A","1.0*@0",RooArgList(n_1_A) );
  RooRealVar   *bfrac = new RooRealVar("bfrac","b fraction", f_init, paramslower[4], paramsupper[4]);
  // Set up crystal ball shapes
  RooCBShape* cb_1_A = new RooCBShape("cball_1_A", "cystal Ball", *(ws->var("mass")), mean, sigma_1_A, alpha_1_A, n_1_A);
  RooAddPdf*  cb_A;
  //DOUBLE CRYSTAL BALL
  RooCBShape* cb_2_A = new RooCBShape("cball_2_A", "cystal Ball", *(ws->var("mass")), mean, sigma_2_A, alpha_2_A, n_2_A);
  cb_A = new RooAddPdf("cb_A","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*bfrac) );
  RooRealVar *nSig= new RooRealVar("nSig","inclusive Jpsi signals",0,100000);
  //RooRealVar *nSigPR= new RooRealVar("nSigPR","prompt Jpsi signals",0,100000);
  //RooRealVar *nSigNP= new RooRealVar("nSigNP","non-prompt Jpsi signals",0,100000);
  cb_A = new RooAddPdf("cb_A","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*bfrac) );

//B:Ctau distribution
  ws->factory("GaussModel::GW_PRRes(Jpsi_Ct,meanPRResW[0.,-0.01,0.01],sigmaPRResW[2.3,1.1,5.5],one[1.0],Jpsi_CtErr)");
  ws->factory("GaussModel::GN_PRRes(Jpsi_Ct,meanPRResN[0.,-0.01,0.01],sigmaPRResN[0.8,0.6,1.2],one,Jpsi_CtErr)");
  ws->factory("AddModel::CtPRRes({GW_PRRes,GN_PRRes},{fracRes[0.01,0.001,0.999]})");

  //RooRealVar dt("dt", "dt", ctauLow, ctauHigh);
  //RooRealVar tau("tau", "tau", 1.548);
  //RooTruthModel tm("tm", "truth model", dt);
  //RooDecay decay_tm("decay_tm", "decay", dt, tau, tm, RooDecay::DoubleSided);//PR Jpsi

  //// Set up total shapes
  ////For the background: SingleSided(Nonprompt), Flipped(??), DoubleSided(Prompt Jpsi)
  //Roo??* ctPR = new Roo??("", "", *(ws->var("ctau3D")), pr);
  //Roo??* ctNP = new Roo??("", "", *(ws->var("ctau3D")), np);
  //Roo??* ctBkg = new Roo??("", "", *(ws->var("ctau3D")), bkg);
  //RooAddPdf* ctTotal;
  //ctTotal = new RooAddPdf("ctTotal","Sig + Bkg ", RooArgList(*ctPR,*ctNP,*ctBkg), RooArgList(*bfrac) );

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

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *bkg = new RooGenericPdf("bkgHighPt_A","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_A));
  RooRealVar *nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",10000,0,5000000); 

  RooFormulaVar fracAlpha("fracAlpha","@0/(@0+@1)",RooArgList(*nBkg,*nSig)); ws->import(fracAlpha);
  //cout<<"Alpha: "<<ws->var("fracAlpha")->getVal()<<endl;

//B: Ctau Bkg function
  ws->factory("Decay::CtBkgPos(Jpsi_Ct,lambdap[0.42,0.1,1.0],CtPRRes,RooDecay::SingleSided)");
  ws->factory("Decay::CtBkgNeg(Jpsi_Ct,lambdam[0.02,0.001,0.1],CtPRRes,RooDecay::Flipped)");
  ws->factory("Decay::CtBkgDbl(Jpsi_Ct,lambdasym[0.2,0.001,0.8],CtPRRes,RooDecay::DoubleSided)");
  ws->factory("SUM::CtBkgSum1(fracCtBkg1[0.9,0.0,1.0]*CtBkgPos,CtBkgNeg)");
  ws->factory("SUM::CtBkgSum2(fracCtBkg2[0.9,0.0,1.0]*CtBkgSum1,CtBkgDbl)");
  ws->factory("SUM::CtBkgTot(fracCtBkg3[0.29,0.0,1.0]*CtPRRes,CtBkgSum2)");

  //Build the model
//Model A: Mass
  RooAddPdf* model_A = new RooAddPdf();
  model_A = new RooAddPdf("model_A","Jpsi + Bkg",RooArgList(*cb_A, *bkg),RooArgList(*nSig,*nBkg));
  //model_A = new RooAddPdf("model_A","PR Jpsi + NP Jpsi + Bkg",RooArgList(*cb_1_A, *cb_2_A, *bkg),RooArgList(*nSigPR,*nSigNP,*nBkg));
  ws->import(*model_A);
  c_A->cd();
  RooPlot* myPlot2_A = (RooPlot*)myPlot_A->Clone();
  ws->data("reducedDS_A")->plotOn(myPlot2_A,Name("dataOS_FIT_A"),MarkerSize(.8));

//Model B: Ctau
  RooAddPdf* model_B = new RooAddPdf();
  model_B = new RooAddPdf("model_B","Jpsi + Bkg",RooArgList(*cb_A, *bkg),RooArgList(*nSig, *nBkg));
  //model_B = new RooAddPdf("model_B","Jpsi + Bkg",RooArgList(*ctPRRes, *CkBkgTot),RooArgList(*nSig, *nBkg));
  ws->import(*model_B);
  c_B->cd();
  RooPlot* myPlot2_B = (RooPlot*)myPlot_B->Clone();
  //ws->data("reducedDS_B")->plotOn(myPlot2_B,Name("dataOS_FIT_B"), Cut(Form("(mass<2.9||mass>3.3)&&%s",accCut.Data())), MarkerSize(.8));
  ws->data("reducedDS_B")->plotOn(myPlot2_B,Name("dataOS_FIT_B"), Cut(Form("mass>%.2f&&mass<%.2f", massLow, massHigh)), MarkerSize(.8));

  // Construct simultaneous PDF
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simPdf",tp);
  simPdf->addPdf(*(ws->pdf("model_A")),"A") ;
  //simPdf->addPdf(*(ws->pdf("model_B")),"B") ;
  ws->import(*simPdf);

  cout << endl << "********* Starting Simutaneous Fit **************" << endl << endl;
  RooFitResult* fitResSim = ws->pdf("simPdf")->fitTo(*dsAB,Save(), Hesse(kTRUE),Range(massLow,massHigh),Timer(kTRUE),Extended(kTRUE));
  //RooFitResult* fitResSim = ws->pdf("simPdf")->fitTo(*dsAB,Save(), Hesse(kTRUE), Timer(kTRUE),Extended(kTRUE));
  cout << endl << "********* Finished Simutaneous Fit **************" << endl << endl;

  c_A->cd();
  gPad->SetLogy();
  ws->pdf("model_A")->plotOn(myPlot2_A,Name("modelHist_A"), LineColor(kBlack));
  ws->pdf("model_A")->plotOn(myPlot2_A,Name("Sig_A"),Components(RooArgSet(*cb_A)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_A")->plotOn(myPlot2_A,Name("bkgPDF_A"),Components(RooArgSet(*bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  myPlot2_A->SetFillStyle(4000);
  myPlot2_A->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2_A->GetYaxis()->SetTitleOffset(1.43);
  //myPlot2_A->GetYaxis()->CenterTitle();
  //myPlot2_A->GetYaxis()->SetTitleSize(0.058);
  myPlot2_A->GetYaxis()->SetLabelSize(0.054);
  //myPlot2_A->GetXaxis()->SetLabelSize(0);
  myPlot2_A->GetXaxis()->SetRangeUser(massLow,massHigh);
  //myPlot2_A->GetYaxis()->SetRangeUser(2*10, 2*10e3);
  myPlot2_A->SetMinimum(2*10);
  //myPlot2_A->GetXaxis()->SetTitleSize(0);
  myPlot2_A->Draw();

  float pos_text_x = 0.15;
  float pos_text_y = 0.816;
  float pos_y_diff = 0.056;
  float text_size = 19;
  int text_color = 1;
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  drawText(Form("n_{J/#psi} = %.f #pm %.f", ws->var("nSig")->getVal(),ws->var("nSig")->getError()),pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  drawText(Form("n_{Bkg} = %.f #pm %.f", ws->var("nBkg")->getVal(),ws->var("nBkg")->getError()),pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
  //drawText(Form("Signal Function : %s CB", SignalCB.Data() ), 0.55,0.54,1,14);

  c_B->cd();
  gPad->SetLogy();
  //ws->pdf("model_B")->plotOn(myPlot2_B,Name("modelHist_B"));
  //ws->pdf("model_B")->plotOn(myPlot2_B,Name("Sig_B"),Components(RooArgSet(*cb_B)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  //ws->pdf("model_B")->plotOn(myPlot2_B,Name("bkgPDF_B"),Components(RooArgSet(*bkg_B)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  myPlot2_B->SetFillStyle(4000);
  //myPlot2_B->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2_B->GetYaxis()->SetTitleOffset(1.43);
  //myPlot2_B->GetYaxis()->CenterTitle();
  //myPlot2_B->GetYaxis()->SetTitleSize(0.058);
  myPlot2_B->GetYaxis()->SetLabelSize(0.054);
  //myPlot2_B->GetXaxis()->SetLabelSize(0);
  myPlot2_B->GetXaxis()->SetRangeUser(ctauLow,ctauHigh);
  myPlot2_B->GetYaxis()->SetRangeUser(10e-2, 10e5);
  //myPlot2_B->GetXaxis()->SetTitleSize(0);
  myPlot2_B->Draw();

  c_A->Update();
  c_B->Update();

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
