#include <iostream>
#include "../../../rootFitHeaders.h"
#include "../../../commonUtility.h"
#include "../../../JpsiUtility.h"
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
#include "../../../cutsAndBin.h"
#include "../../../CMS_lumi_v2mass.C"
#include "../../../tdrstyle.C"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit;
void CtauErr_v2Bin_n(
    double v2=0,
    float ptLow=3, float ptHigh=4.5,
    float yLow=1.6, float yHigh=2.4,
    int cLow=20, int cHigh=120,
    int PRw=1, bool fEffW = true, bool fAccW = true, bool isPtW = true, bool isTnP = true
    )
{

  //TString DATE="210507";
  //TString DATE="10_60";
  //TString DATE="20_40";
  //TString DATE="0_180";
  TString DATE;
  if(ptLow==6.5&&ptHigh==50&&!(cLow==0&&cHigh==180)) DATE=Form("%i_%i",0,180);
  else DATE=Form("%i_%i",cLow/2,cHigh/2);
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir(Form("roots/2DFit_%s/CtauErr_v2Bins_n",DATE.Data()),kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/CtauErr_v2Bins_n",DATE.Data()),kTRUE);

  TString fname;
  if (PRw==1) fname="PR";
  else if (PRw==2) fname="NP";

  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  TFile* f1; TFile* f2; TFile* f3; TFile* fMass;
  TString kineCut;
  TString v2Cut;
  TString SigCut;
  TString BkgCut;
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);

  f1 = new TFile(Form("../../../skimmedFiles/v2Cut_Nom/OniaRooDataSet_isMC0_JPsi_%s_m2.6-3.5_OS_Effw%d_Accw%d_PtW%d_TnP%d_%s_220113.root",kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,DATE.Data()));
  fMass = new TFile(Form("roots/2DFit_%s/Mass_v2Bins_n/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_v2_%.2f.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP, v2));
  kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.5 && cBin>%d && cBin<%d",ptLow, ptHigh, yLow, yHigh, cLow, cHigh);
  if (v2==2.1) v2Cut = Form("&& v2>%.2f && v2<%.2f",v2,v2+0.6);
  else if (ptLow==25&&ptHigh==50&&v2==1.8)  v2Cut = Form("&& v2>%.2f && v2<%.2f",v2,v2+0.9);
  else if (ptLow==25&&ptHigh==50&&v2==-2.7)  v2Cut = Form("&& v2>%.2f && v2<%.2f",v2,v2+0.9);
  else if (ptLow==25&&ptHigh==50&&v2==-1.8)  v2Cut = Form("&& v2>%.2f && v2<%.2f",v2,v2+0.6);
  else if (ptLow==25&&ptHigh==50&&v2==1.2)  v2Cut = Form("&& v2>%.2f && v2<%.2f",v2,v2+0.6);
  else if (ptLow==20&&ptHigh==25&&v2==1.8)  v2Cut = Form("&& v2>%.2f && v2<%.2f",v2,v2+0.9);
  else if (ptLow==20&&ptHigh==25&&v2==-2.7)  v2Cut = Form("&& v2>%.2f && v2<%.2f",v2,v2+0.9);
  else if (v2==-2.7) v2Cut = Form("&& v2>%.2f && v2<%.2f",v2,v2+0.6);
  else if (v2==-10) v2Cut = Form("&& v2>%.2f && v2<%.2f",v2,v2+20.00);
  else v2Cut = Form("&& v2>%.2f && v2<%.2f", v2, v2+0.3);

  cout << "v2 Cut : " << v2Cut << endl;

  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

  TString OS="recoQQsign==0 &&";

  kineCut = OS+accCut+kineCut+v2Cut;

  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooDataSet *datasetMass = (RooDataSet*)fMass->Get("datasetMass");
  RooAddPdf* pdfMASS_Tot = (RooAddPdf*)fMass->Get("pdfMASS_Tot");

  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  ws->import(*datasetMass);
  ws->import(*pdfMASS_Tot);
  RooArgSet *argSet = new RooArgSet( *(ws->var("ctau3D")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("weight")), *(ws->var("ctau3DRes")), *(ws->var("ctau3DErr")) );
  argSet->add(*(ws->var("pt1")) ); argSet->add(*(ws->var("pt2")) ); argSet->add(*(ws->var("eta1")) );  argSet->add(*(ws->var("eta2")) ); argSet->add(*(ws->var("recoQQsign")) ); argSet->add(*(ws->var("cBin"))); argSet->add(*(ws->var("v2")));
  RooDataSet *datasetW = new RooDataSet("datasetW","A sample",
		  *argSet,
		  Import(*dataset),WeightVar(*ws->var("weight")));
  RooDataSet *datasetWo = new RooDataSet("datasetWo","A sample",
		  *argSet,
		  Import(*dataset));
  ws->import(*datasetW);
  ws->import(*datasetWo);

  RooDataSet *dsAB = (RooDataSet*)datasetW->reduce(*argSet, kineCut.Data() );
  RooDataSet *dsAB1 = (RooDataSet*)datasetWo->reduce(*argSet, kineCut.Data() );
  dsAB->Print();
  cout << "Weight : " << ws->var("weight")->getVal() << endl;
  cout << "pt: "<<ptLow<<"-"<<ptHigh<<", y: "<<yLow<<"-"<<yHigh<<", Cent: "<<cLow<<"-"<<cHigh<<"%"<<endl;
  cout << "####################################" << endl;
  dsAB->SetName("dsAB");
  dsAB1->SetName("dsAB1");
  ws->import(*dsAB);
  ws->import(*dsAB1);

  //ws->var("ctau3DErr")->setRange(0,0.25);
    ws->var("ctau3DErr")->setRange(ctauErrLow, ctauHigh);
  /*
  TCanvas* c_D =  new TCanvas("canvas_D","My plots",4,420,550,520);
  c_D->cd();
  TPad *pad_D_1 = new TPad("pad_D_1", "pad_D_1", 0, 0.16, 0.98, 1.0);
  pad_D_1->SetTicks(1,1);
  pad_D_1->Draw(); pad_D_1->cd();

  RooPlot* myPlot_D = ws->var("ctau3DErr")->frame(Bins(80), Range(0, 0.25));
  RooPlot* myPlot2_D = (RooPlot*)myPlot_D->Clone();
//  ws->data("datasetW")->plotOn(myPlot2_D,Name("MCHist_Tot_NoW"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlue+2), MarkerSize(.7));
//  ws->data("datasetWo")->plotOn(myPlot2_D,Name("MCHist_Tot"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7));
  ws->data("dsAB")->plotOn(myPlot2_D,Name("MCHist_Tot1"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), MarkerColor(kRed+2));
  ws->data("dsAB1")->plotOn(myPlot2_D,Name("MCHist_Tot2"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), MarkerColor(kBlue+2));
  myPlot2_D->Draw();*/
  //int nBins = min(int( round((ctauErrHigh - ctauErrLow)/0.0025) ), 100);
  int nBins = min(int( round((ws->var("ctau3DErr")->getMax() - ws->var("ctau3DErr")->getMin())/0.0025) ), 100);
    cout << " nBin : " << nBins << endl;
//  RooBinning bins(nBins, ctauErrLow, ctauErrHigh);
  //RooRealVar *N_Jpsi = ws->var("N_Jpsi");
  //RooRealVar *N_Bkg = ws->var("N_Bkg");
  //ws->var("ctau3DErr")->setRange(ctauErrLow, ctauErrHigh);
  //ws->var("ctau3DErr")->setRange(ctauErrLow,ctauErrHigh,ctauErrLow, ctauErrHigh);
  //ws->var("ctau3DErr")->Print();
  //ws->var("N_Jpsi")->Print();
  //***********************************************************************
  //**************************** CTAU ERR FIT *****************************
  //***********************************************************************
  cout << endl << "*************** Start SPLOT *****************" << endl << endl;
  //SPlot Ctau Error
  RooRealVar *sigYield = ws->var("N_Jpsi");
  RooRealVar *bkgYield = ws->var("N_Bkg");
  RooArgList yieldList;
  yieldList.add(*ws->var("N_Jpsi"));
  yieldList.add(*ws->var("N_Bkg"));
  cout<<"Sig Yield: "<<sigYield->getVal()<<" +/- "<<sigYield->getError()<<endl;
  cout<<"Bkg Yield: "<<bkgYield->getVal()<<" +/- "<<bkgYield->getError()<<endl;
  RooDataSet* data = (RooDataSet*)ws->data("dsAB")->Clone("TMP_DATA");
  RooArgSet* cloneSet = (RooArgSet*)RooArgSet(*ws->pdf("pdfMASS_Tot"),"pdfMASS_Tot").snapshot(kTRUE);
  auto clone_mass_pdf = (RooAbsPdf*)cloneSet->find("pdfMASS_Tot");
  clone_mass_pdf->setOperMode(RooAbsArg::ADirty, kTRUE);

  RooStats::SPlot sData = RooStats::SPlot("sData", "An SPlot", *data, clone_mass_pdf, yieldList);
  ws->import(*data, Rename("dataset_SPLOT"));
  cout<<"[INFO] Jpsi yield -> Mass Fit:"<<ws->var("N_Jpsi")->getVal()<<", sWeights :"<<sData.GetYieldFromSWeight("N_Jpsi")<<endl;
  cout<<"[INFO] Bkg  yield -> Mass Fit:"<<ws->var("N_Bkg")->getVal()<<", sWeights :"<<sData.GetYieldFromSWeight("N_Bkg")<<endl;
  //create weighted data sets
  //total
  RooPlot* myPlot_B = ws->var("ctau3DErr")->frame(Range(ctauErrLow,ctauErrHigh));
  TH1D* hTot = (TH1D*)ws->data("dsAB")->createHistogram(("hTot"), *ws->var("ctau3DErr"),
		  //Binning(myPlot_B->GetNbinsX(), myPlot_B->GetXaxis()->GetXmin(), myPlot_B->GetXaxis()->GetXmax()));
          Binning(myPlot_B->GetNbinsX(), myPlot_B->GetXaxis()->GetXmin(), myPlot_B->GetXaxis()->GetXmax()));
  double ctauErrMin; double ctauErrMax;
  TH1D* hTot_M = (TH1D*)ws->data("dsAB")->createHistogram(("hTot_M"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrLow,ctauErrHigh));
  //Bkg
  RooDataSet* dataw_Bkg_b = new RooDataSet("dataw_Bkg_b","TMP_BKG_DATA", (RooDataSet*)ws->data("dataset_SPLOT"),
      RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Bkg_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Bkg_sw");
  TH1D* hBkg = (TH1D*)dataw_Bkg_b->createHistogram(("hBkg"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrLow,ctauErrHigh));

  //data
  RooDataSet* dataw_Sig_b = new RooDataSet("dataw_Sig_b","TMP_SIG_DATA", (RooDataSet*)ws->data("dataset_SPLOT"),
      RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Jpsi_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Jpsi_sw");
  TH1D* hSig = (TH1D*)dataw_Sig_b->createHistogram(("hSig"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrLow,ctauErrHigh));

  for(int i=0; i<hTot_M->GetNbinsX()/2; i++){
    //if(hSig->GetBinContent(i)<=0&&hSig->GetBinContent(i+1)<=0&&hSig->GetBinContent(i+2)>=1&&hSig->GetBinContent(i+3)>=1){
    if(hTot_M->GetBinContent(i)>1){//pt 7-7.5
      //if(ptLow==3&&ptHigh==4.5&&cLow==20&&cHigh==120) ctauErrMin = hSig->GetBinLowEdge(i);
      //else ctauErrMin = hSig->GetBinLowEdge(i+1);
      ctauErrMin = hTot_M->GetBinLowEdge(i+1);
      break;}
  }
  
  //if(ptLow==3&&ptHigh==4.5&&cLow==20&&cHigh==120) ctauErrMin = 0;

  //for(int i=hSig->GetNbinsX()/2; i<hSig->GetNbinsX(); i++){
  //for(int i=20; i<hSig->GetNbinsX(); i++){
  //  if(hSig->GetBinContent(i)>=1&&hSig->GetBinContent(i+1)<1&&hSig->GetBinContent(i+2)<1){
  for(int i=0; i<hTot_M->GetNbinsX(); i++){
    if(hTot_M->GetBinContent(i)>=1&&hTot_M->GetBinContent(i+1)<1&&hTot_M->GetBinContent(i+2)<1){
      //ctauErrMax = hSig->GetBinLowEdge(i)+hSig->GetBinWidth(i); 
      ctauErrMax = hTot_M->GetBinLowEdge(i)+hTot_M->GetBinWidth(i);
      break;}
    //else { ctauErrMax = hSig->GetBinLowEdge(i); }
    else { ctauErrMax = hTot_M->GetBinLowEdge(i); }
  }

  if(v2==-0.9&&ptLow==3&&ptHigh==4.5) ctauErrMax=0.1566;
  else if(v2==-10) ctauErrMax=0.1293;
  else if(v2==-0.6&&ptLow==3&&ptHigh==4.5) ctauErrMax=0.1404;
  else if(v2==-1.8&&ptLow==3&&ptHigh==4.5) ctauErrMax=0.113;
  else if(v2==-2.1&&ptLow==3&&ptHigh==4.5) ctauErrMax=0.0954;
  else if(v2==-2.4&&ptLow==3&&ptHigh==4.5) ctauErrMax=0.099;
  else if(v2==0.6&&ptLow==3&&ptHigh==4.5) ctauErrMax=0.1547;
  else if(v2==0.9&&ptLow==3&&ptHigh==4.5) ctauErrMax=0.1455;
  else if(v2==1.2&&ptLow==3&&ptHigh==4.5) ctauErrMax=0.1074;
  else if(v2==1.5&&ptLow==3&&ptHigh==4.5) ctauErrMax=0.1277;
  else if(v2==1.8&&ptLow==3&&ptHigh==4.5) ctauErrMax=0.1152;
  else if(v2==2.1&&ptLow==3&&ptHigh==4.5) ctauErrMax=0.0932;
  else if(v2==2.4&&ptLow==3&&ptHigh==4.5) ctauErrMax=0.0882;

  else if(v2==-0.3&&ptLow==3&&ptHigh==6.5) ctauErrMax=0.1363;
  else if(v2==-0.9&&ptLow==3&&ptHigh==6.5) ctauErrMax=0.1566;
  else if(v2==-2.7&&ptLow==3&&ptHigh==6.5) ctauErrMax=0.0878;
  else if(v2==0.3&&ptLow==3&&ptHigh==6.5) ctauErrMax=0.1656;
  else if(v2==0.6&&ptLow==3&&ptHigh==6.5) ctauErrMax=0.1476;
  else if(v2==0.9&&ptLow==3&&ptHigh==6.5) ctauErrMax=0.1382;
  else if(v2==1.5&&ptLow==3&&ptHigh==6.5) ctauErrMax=0.1152;
  else if(v2==1.8&&ptLow==3&&ptHigh==6.5) ctauErrMax=0.1098;

  else if(v2==-0.3&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.1547;
  else if(v2==-0.6&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.133;
  else if(v2==-0.9&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.1296;
  else if(v2==-1.2&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.126;
  else if(v2==-1.5&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.1224;
  else if(v2==-1.8&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.1224;
  else if(v2==-2.1&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.0972;
  else if(v2==0&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.1149;
  else if(v2==0.3&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.1237;
  else if(v2==0.6&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.1242;
  else if(v2==0.9&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.1314;
  else if(v2==1.5&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.1224;
  else if(v2==1.8&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.0824;
  else if(v2==2.1&&ptLow==4.5&&ptHigh==6.5) ctauErrMax=0.0648;

  else if(v2==-0.3&&ptLow==6.5&&ptHigh==7.5) ctauErrMax=0.0642;
  else if(v2==0.6&&ptLow==6.5&&ptHigh==7.5) ctauErrMax=0.1152;
  else if(v2==1.5&&ptLow==6.5&&ptHigh==7.5) ctauErrMax=0.0665;
  else if(v2==-2.1&&ptLow==6.5&&ptHigh==7.5) ctauErrMax=0.0646;
  else if(v2==-2.7&&ptLow==6.5&&ptHigh==7.5) ctauErrMax=0.0626;
  else if(v2==2.1&&ptLow==6.5&&ptHigh==7.5) ctauErrMax=0.0558;

  else if(v2==-1.5&&ptLow==6.5&&ptHigh==7.5) ctauErrMax=0.1062;

  else if(v2==-1.5&&ptLow==6.5&&ptHigh==9) ctauErrMax=0.1026;
  else if(v2==-1.8&&ptLow==6.5&&ptHigh==9) ctauErrMax=0.0916;
  else if(v2==-2.7&&ptLow==6.5&&ptHigh==9) ctauErrMax=0.0474;
  else if(v2==-0.3&&ptLow==6.5&&ptHigh==9) ctauErrMax=0.102;
  else if(v2==0.3&&ptLow==6.5&&ptHigh==9) ctauErrMax=0.1314;
  else if(v2==0.6&&ptLow==6.5&&ptHigh==9) ctauErrMax=0.1059;
  else if(v2==1.2&&ptLow==6.5&&ptHigh==9) ctauErrMax=0.106;
  else if(v2==2.1&&ptLow==6.5&&ptHigh==9) ctauErrMax=0.0648;

  else if(v2==-0.3&&ptLow==7.5&&ptHigh==9) ctauErrMax=0.099;
  else if(v2==-1.5&&ptLow==7.5&&ptHigh==9) ctauErrMax=0.0738;
  else if(v2==0&&ptLow==7.5&&ptHigh==9) ctauErrMax=0.1044;
  else if(v2==0.6&&ptLow==7.5&&ptHigh==9) ctauErrMax=0.0918;
  else if(v2==0.9&&ptLow==7.5&&ptHigh==9) ctauErrMax=0.0933;

  else if(v2==-1.2&&ptLow==9&&ptHigh==12) ctauErrMax=0.0702;
  else if(v2==-1.8&&ptLow==9&&ptHigh==12) ctauErrMax=0.0648;
  else if(v2==2.1&&ptLow==9&&ptHigh==12) ctauErrMax=0.0396;
  
  else if(v2==-0.6&&ptLow==12&&ptHigh==15) ctauErrMax=0.0518;
  else if(v2==-0.9&&ptLow==12&&ptHigh==15) ctauErrMax=0.0648;

  else if(v2==0&&ptLow==15&&ptHigh==20) ctauErrMax=0.0468;
  else if(v2==0.9&&ptLow==15&&ptHigh==20) ctauErrMax=0.054;
  else if(v2==-2.7&&ptLow==15&&ptHigh==20) ctauErrMax=0.027;

  else if(v2==-1.5&&ptLow==15&&ptHigh==50) ctauErrMax=0.0398;
  else if(v2==0.3&&ptLow==15&&ptHigh==50) ctauErrMax=0.0504;

  else if(v2==0.3&&ptLow==20&&ptHigh==25) ctauErrMax=0.0342;
  else if(v2==1.8&&ptLow==20&&ptHigh==25) ctauErrMax=0.0234;
  else if(v2==-2.7&&ptLow==20&&ptHigh==25) ctauErrMax=0.0234;

  else if(v2==0.3&&ptLow==20&&ptHigh==50) ctauErrMax=0.0342;

  else if(v2==0&&ptLow==25&&ptHigh==50) ctauErrMax=0.0306;

  else if(v2==-0.3&&cLow==0&&cHigh==20) ctauErrMax=0.0915;
  else if(v2==-0.6&&cLow==0&&cHigh==20) ctauErrMax=0.1008;
  else if(v2==-0.9&&cLow==0&&cHigh==20) ctauErrMax=0.1062;
  else if(v2==-1.2&&cLow==0&&cHigh==20) ctauErrMax=0.1026;
  else if(v2==-1.5&&cLow==0&&cHigh==20) ctauErrMax=0.0828;
  else if(v2==-2.7&&cLow==0&&cHigh==20) { ctauErrMin=0.0126; ctauErrMax=0.0684; }
  else if(v2==0&&cLow==0&&cHigh==20) ctauErrMax=0.1134;
  else if(v2==0.3&&cLow==0&&cHigh==20) ctauErrMax=0.09;
  else if(v2==0.6&&cLow==0&&cHigh==20) ctauErrMax=0.1044;
  else if(v2==1.5&&cLow==0&&cHigh==20) ctauErrMax=0.07536;

  else if(v2==-0.3&&cLow==20&&cHigh==40) ctauErrMax=0.0894;
  else if(v2==-1.2&&cLow==20&&cHigh==40) ctauErrMax=0.0972;
  else if(v2==-2.1&&cLow==20&&cHigh==40) ctauErrMax=0.0666;
//  else if(v2==-2.7&&cLow==20&&cHigh==40) ctauErrMax=0.0468;
  else if(v2==-1.5&&cLow==20&&cHigh==40) ctauErrMax=0.0808;
  else if(v2==0.6&&cLow==20&&cHigh==40) ctauErrMax=0.1062;
  else if(v2==1.5&&cLow==20&&cHigh==40) ctauErrMax=0.0794;

  else if(v2==-1.2&&cLow==40&&cHigh==60) ctauErrMax=0.0774;
  else if(v2==-0.6&&cLow==40&&cHigh==60) ctauErrMax=0.0952;
  else if(v2==0.3&&cLow==40&&cHigh==60) ctauErrMax=0.108;
  else if(v2==0.6&&cLow==40&&cHigh==60) ctauErrMax=0.0918;
  else if(v2==0.9&&cLow==40&&cHigh==60) ctauErrMax=0.0972;
  else if(v2==1.8&&cLow==40&&cHigh==60) ctauErrMax=0.0646;

  else if(v2==-0.3&&cLow==60&&cHigh==80) ctauErrMax=0.0864;
  else if(v2==-0.6&&cLow==60&&cHigh==80) ctauErrMax=0.097;
  else if(v2==-1.2&&cLow==60&&cHigh==80) ctauErrMax=0.0756;
  else if(v2==-1.5&&cLow==60&&cHigh==80) ctauErrMax=0.072;
  else if(v2==-1.8&&cLow==60&&cHigh==80) ctauErrMax=0.0498;
  else if(v2==-2.1&&cLow==60&&cHigh==80) ctauErrMin=0.0144;
  else if(v2==0.6&&cLow==60&&cHigh==80) ctauErrMax=0.0972;
  else if(v2==0.9&&cLow==60&&cHigh==80) ctauErrMax=0.0846;

  else if(v2==-0.3&&cLow==80&&cHigh==100) ctauErrMax=0.0648;
  else if(v2==-0.6&&cLow==80&&cHigh==100) ctauErrMax=0.0666;
  else if(v2==-0.9&&cLow==80&&cHigh==100) ctauErrMax=0.0648;
  else if(v2==-1.2&&cLow==80&&cHigh==100) ctauErrMax=0.0738;
  else if(v2==-1.8&&cLow==80&&cHigh==100) ctauErrMax=0.045;
  else if(v2==2.1&&cLow==80&&cHigh==100) ctauErrMax=0.0504;
  else if(v2==0&&cLow==80&&cHigh==100) ctauErrMax=0.0882;
  else if(v2==0.3&&cLow==80&&cHigh==100) ctauErrMax=0.0603;
  else if(v2==1.5&&cLow==80&&cHigh==100) ctauErrMax=0.0648;

  else if(v2==-1.8&&cLow==100&&cHigh==180) ctauErrMax=0.0666;
  else if(v2==-1.5&&cLow==100&&cHigh==180) ctauErrMax=0.063;
  else if(v2==-2.7&&cLow==100&&cHigh==180) ctauErrMax=0.054;
  else if(v2==0&&cLow==100&&cHigh==180) ctauErrMax=0.0828;
  else if(v2==1.8&&cLow==100&&cHigh==180) { ctauErrMin=0.01184; ctauErrMax=0.0523; }
 /* 
  else if(v2==-0.6&&ptLow==6.5&&ptHigh==7.5) ctauErrMax=0.0918;
  else if(v2==-0.9&&ptLow==6.5&&ptHigh==7.5) ctauErrMax=0.0792;
  else if(v2==0&&ptLow==6.5&&ptHigh==7.5) ctauErrMax=0.0684;
  else if(v2==0.6&&ptLow==6.5&&ptHigh==7.5) ctauErrMax=0.072;

  else if(v2==-0.6&&ptLow==7.5&&ptHigh==9) ctauErrMax=0.0792;
  else if(v2==-2.7&&ptLow==7.5&&ptHigh==9) ctauErrMax=0.0846;
  else if(v2==0&&ptLow==7.5&&ptHigh==9) ctauErrMax=0.0846;
  else if(v2==0.6&&ptLow==7.5&&ptHigh==9) ctauErrMax=0.081; //0.0774
  else if(v2==1.5&&ptLow==7.5&&ptHigh==9) ctauErrMax=0.0774; //0.0702
  else if(v2==2.4&&ptLow==7.5&&ptHigh==9) ctauErrMax=0.0756; //0.072

  else if(v2==-2.1&&ptLow==15&&ptHigh==50) ctauErrMax=0.0486; //0.072
*/
  //if(ptLow==12&&ptHigh==15&&cLow==20&&cHigh==120) ctauErrMax=0.05;
  //else if(ptLow==6.5&&ptHigh==50&&cLow==0&&cHigh==20) ctauErrMax=0.1;
  //else if(ptLow==6.5&&ptHigh==50&&cLow==20&&cHigh==40) ctauErrMax=0.105;
  //else if(ptLow==6.5&&ptHigh==50&&cLow==40&&cHigh==60) ctauErrMax=0.1;
  //else if(ptLow==6.5&&ptHigh==50&&cLow==60&&cHigh==80) ctauErrMax=0.1;
  //else if(ptLow==6.5&&ptHigh==50&&cLow==80&&cHigh==100) ctauErrMax=0.09;
  //else if(ptLow==6.5&&ptHigh==50&&cLow==100&&cHigh==180) ctauErrMax=0.065;

  cout << "ctauErrMax : " << ctauErrMax << " ctauErrMin : " << ctauErrMin << endl;

    double BinWidth = (ctauErrHigh-ctauErrLow)/nBins;
    cout << " BinWidth : " << BinWidth << endl;
    int newBins = (ctauErrMax - ctauErrMin)/BinWidth ;
    cout << " newBins : " << newBins << endl;
  //TH1D* hTot_w = (TH1D*)ws->data("dsAB")->createHistogram(("hTot_w"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrMin,ctauErrMax));

  //cout<<"\n Search Min"<<endl;
  //for(int i=0; i<hTot_w->GetNbinsX()/2; i++){
  //    cout<<"Bin: "<<i<<", CtauError: "<<hTot_w->GetBinLowEdge(i)<<", Events: "<<hTot_w->GetBinContent(i)<<endl;
  //  if(hTot_w->GetBinContent(i)<=0&&hTot_w->GetBinContent(i+1)<=0&&hTot_w->GetBinContent(i+2)>=1&&hTot_w->GetBinContent(i+3)>=1){
  //    cout<<"##### i: "<<i+1<<", CtauError:  "<<hTot_w->GetBinLowEdge(i+1)<<hTot_w->GetBinContent(i+1)<<", Events: "<<hTot_w->GetBinContent(i+1)<<endl;
  //    cout<<"##### i(Min): "<<i+2<<", CtauError:  "<<hTot_w->GetBinLowEdge(i+2)<<hTot_w->GetBinContent(i+2)<<", Events: "<<hTot_w->GetBinContent(i+2)<<endl;
  //    cout<<"##### i: "<<i+3<<", CtauError:  "<<hTot_w->GetBinLowEdge(i+3)<<hTot_w->GetBinContent(i+3)<<", Events: "<<hTot_w->GetBinContent(i+3)<<endl;
  //    ctauErrMin = hTot_w->GetBinLowEdge(i+2);
  //    break;}
  //}

  //cout<<"\n Search Max"<<endl;
  //for(int i=20; i<hTot_w->GetNbinsX(); i++){
  //    cout<<"Bin: "<<i<<", CtauError: "<<hTot_w->GetBinLowEdge(i)+hTot_w->GetBinWidth(i)<<", Events: "<<hTot_w->GetBinContent(i)<<endl;
  //  if(hTot_w->GetBinContent(i)>=1&&hTot_w->GetBinContent(i+1)<1){
  //  //if(hTot_w->GetBinContent(i+1)<1&&hTot_w->GetBinContent(i+2)<1){
  //  //if(hTot_w->GetBinContent(i)>=1&&hTot_w->GetBinContent(i+1)<1&&hTot_w->GetBinContent(i+2)<1){
  //  //if(hTot_w->GetBinContent(i)>=100&&hTot_w->GetBinContent(i+1)<50&&hTot_w->GetBinContent(i+2)<50){
  //    cout<<"##### i(Max): "<<i<<", CtauError:  "<<hTot_w->GetBinLowEdge(i)+hTot_w->GetBinWidth(i)<<hTot_w->GetBinContent(i)<<", Events: "<<hTot_w->GetBinContent(i)<<endl;
  //    cout<<"##### i: "<<i+1<<", CtauError:  "<<hTot_w->GetBinLowEdge(i+1)+hTot_w->GetBinWidth(i+1)<<hTot_w->GetBinContent(i+1)<<", Events: "<<hTot_w->GetBinContent(i+1)<<endl;
  //    cout<<"##### i: "<<i+2<<", CtauError:  "<<hTot_w->GetBinLowEdge(i+2)+hTot_w->GetBinWidth(i+2)<<hTot_w->GetBinContent(i+2)<<", Events: "<<hTot_w->GetBinContent(i+2)<<endl;
  //    ctauErrMax = hTot_w->GetBinLowEdge(i-1)+hTot_w->GetBinWidth(i-1);
  //    break;}
  //  else { ctauErrMax = hTot_M->GetBinLowEdge(i); }
  //}

  
  //TH1D* test = (TH1D*)ws->data("dsAB")->createHistogram(("test"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrMin,ctauErrMax));
  TH1D* test = (TH1D*)ws->data("dsAB")->createHistogram(("test"), *ws->var("ctau3DErr"),Binning(newBins,ctauErrMin,ctauErrMax));
  RooDataHist* totHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), test);
  RooHistPdf* pdfCTAUERR_Tot = new RooHistPdf("pdfCTAUERR_Tot","hist pdf", *ws->var("ctau3DErr"), *totHist);
  //hTot_w->Draw();

  //Bkg
  //TH1D* hBkg_w = (TH1D*)dataw_Bkg_b->createHistogram(("hBkg_w"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrMin,ctauErrMax));
  TH1D* hBkg_w = (TH1D*)dataw_Bkg_b->createHistogram(("hBkg_w"), *ws->var("ctau3DErr"),Binning(newBins,ctauErrMin,ctauErrMax));
  RooDataHist* bkgHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hBkg_w);
  RooHistPdf* pdfCTAUERR_Bkg = new RooHistPdf("pdfCTAUERR_Bkg","hist pdf", *ws->var("ctau3DErr"), *bkgHist);
  //Data
  //TH1D* hSig_w = (TH1D*)dataw_Sig_b->createHistogram(("hSig_w"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrMin,ctauErrMax));
  TH1D* hSig_w = (TH1D*)dataw_Sig_b->createHistogram(("hSig_w"), *ws->var("ctau3DErr"),Binning(newBins,ctauErrMin,ctauErrMax));
  RooDataHist* sigHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hSig_w);
  RooHistPdf* pdfCTAUERR_Jpsi = new RooHistPdf("pdfCTAUERR_Jpsi","hist pdf", *ws->var("ctau3DErr"), *sigHist);
  //Re dataset
  RooDataSet* dataw_Bkg = new RooDataSet("dataw_Bkg","TMP_BKG_DATA", (RooDataSet*)ws->data("dataset_SPLOT"),
      RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Bkg_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass"), *ws->var("v2")), 0, "N_Bkg_sw");
  RooDataSet* dataw_Sig = new RooDataSet("dataw_Sig","TMP_SIG_DATA", (RooDataSet*)ws->data("dataset_SPLOT"),
      RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Jpsi_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass"), *ws->var("v2")), 0, "N_Jpsi_sw");

  //import
  ws->import(*dataw_Sig);
  ws->import(*dataw_Bkg);
  ws->import(*pdfCTAUERR_Tot);
  ws->import(*pdfCTAUERR_Jpsi);
  ws->import(*pdfCTAUERR_Bkg);

  double minRange = (double)(floor(ctauErrMin*100.)/100.);
  double maxRange = (double)(ceil(ctauErrMax*100.)/100.);
  ws->var("ctau3DErr")->setRange("ctauErrWindow", ctauErrMin, ctauErrMax);
  myPlot_B = ws->var("ctau3DErr")->frame(Bins(nBins),Range(minRange-0.01, maxRange+0.01)); // modified

  cout<<ws->data("dsAB")->numEntries()<<endl;
  cout<<ws->data("dataset_SPLOT")->numEntries()<<endl;

  TCanvas* c_B =  new TCanvas("canvas_B","My plots",554,4,550,520);
  c_B->cd();
  TPad *pad_B_1 = new TPad("pad_B_1", "pad_B_1", 0, 0.16, 0.98, 1.0);
  pad_B_1->SetTicks(1,1);
  pad_B_1->Draw(); pad_B_1->cd();
  myPlot_B->SetTitle("");

  c_B->cd();
  c_B->SetLogy();

  pad_B_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot2_B = (RooPlot*)myPlot_B->Clone();
  //ws->data("dataset_SPLOT")->plotOn(myPlot2_B,Name("dataCTAUERR_Tot"), MarkerSize(.7), Binning(nBins));//Normalization(wsmc->data("reducedDS_MC")->sumEntries()
  //ws->data("dsAB")->plotOn(myPlot2_B,Name("dataCTAUERR_Tot"), MarkerSize(.7), Binning(nBins));//Normalization(wsmc->data("reducedDS_MC")->sumEntries()
  ws->data("dsAB")->plotOn(myPlot2_B,Name("dataCTAUERR_Tot"), MarkerSize(.7), Binning(newBins));//Normalization(wsmc->data("reducedDS_MC")->sumEntries()
  ws->pdf("pdfCTAUERR_Tot")->plotOn(myPlot2_B,Name("pdfCTAUERR_Tot"), LineColor(kGreen+1), Range(ctauErrMin,ctauErrMax), LineWidth(2),Normalization(ws->data("dsAB")->sumEntries(),RooAbsReal::NumEvent));
  //ws->data("dataw_Sig")->plotOn(myPlot2_B,Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2), Binning(nBins));
  ws->data("dataw_Sig")->plotOn(myPlot2_B,Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2), Binning(newBins));
  ws->pdf("pdfCTAUERR_Jpsi")->plotOn(myPlot2_B,Name("pdfCTAUERR_Jpsi"),LineColor(kRed+2), LineWidth(2), Range(ctauErrMin,ctauErrMax));
  //ws->data("dataw_Bkg")->plotOn(myPlot2_B,Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue+2), MarkerColor(kBlue+2), Binning(nBins));
  ws->data("dataw_Bkg")->plotOn(myPlot2_B,Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue+2), MarkerColor(kBlue+2), Binning(newBins));
  ws->pdf("pdfCTAUERR_Bkg")->plotOn(myPlot2_B,Name("pdfCTAUERR_Bkg"), LineColor(kBlue+2), LineWidth(2), Range(ctauErrMin,ctauErrMax));
  Double_t YMax = hTot->GetBinContent(hTot->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i=1; i<=hTot->GetNbinsX(); i++) if (hTot->GetBinContent(i)>0) YMin = min(YMin, hTot->GetBinContent(i));
  Double_t Yup(0.),Ydown(0.);
  Yup = YMax*TMath::Power((YMax/YMin), (0.4/(1.0-0.4-0.3)));
  Ydown = YMin/(TMath::Power((YMax/YMin), (0.3/(1.0-0.4-0.3))));
  myPlot2_B->GetYaxis()->SetRangeUser(Ydown,Yup);
  //cout<<ctauErrLow<<", "<<ctauErrHigh<<endl;
  cout<<ws->var("ctau3DErr")->getMin()<<", "<<ws->var("ctau3DErr")->getMax()<<endl;
    RooDataSet *ctauResCutDS =(RooDataSet*)dataw_Sig->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")),*(ws->var("ctau3DErr"))),Form("ctau3DErr>%f&&ctau3DErr<%f",ctauErrMin,ctauErrMax));
    ctauResCutDS->SetName("ctauResCutDS");
    ws->import(*ctauResCutDS);
    Double_t outTot = ws->data("dsAB")->numEntries();
    Double_t outRes = ws->data("ctauResCutDS")->numEntries();
    cout<<"Tot evt: ("<<outTot<<")"<<endl;
    cout<<"Res evt: ("<<outRes<<")"<<endl;
    cout<<"lost evt: ("<<((outTot-outRes)*100)/outTot<<")%, "<<outRes<<"evts"<<endl;

  TLine   *minline = new TLine(ctauErrMin, 0.0, ctauErrMin, (Ydown*TMath::Power((Yup/Ydown),0.5)));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  myPlot2_B->addObject(minline);
  TLine   *maxline = new TLine(ctauErrMax, 0.0, ctauErrMax, (Ydown*TMath::Power((Yup/Ydown),0.5)));
  maxline->SetLineStyle(2);
  maxline->SetLineColor(1);
  maxline->SetLineWidth(3);
  myPlot2_B->addObject(maxline);
  
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
  leg_B->AddEntry(myPlot2_B->findObject("dataCTAUERR_Tot"),"Data","pe");
  leg_B->AddEntry(myPlot2_B->findObject("pdfCTAUERR_Tot"),"Total PDF","l");
  leg_B->AddEntry(myPlot2_B->findObject("pdfCTAUERR_Jpsi"),"Signal","l");
  leg_B->AddEntry(myPlot2_B->findObject("pdfCTAUERR_Bkg"),"Background","l");
  leg_B->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);
  //drawText(Form("n_{J/#psi} = %.f #pm %.f",sData.GetYieldFromSWeight("N_Jpsi"),ws->var("N_Jpsi")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
  //drawText(Form("n_{Bkg} = %.f #pm %.f",sData.GetYieldFromSWeight("N_Bkg"), ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
  drawText(Form("Loss: (%.4f%s) %.f evts", (outTot-outRes)*100/outTot, "%", outTot-outRes),text_x,text_y-y_diff*3,text_color,text_size);
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
  RooHist* hpull_B = frameTMP_B->pullHist("dataCTAUERR_Tot","pdfCTAUERR_Tot");
  //RooHist* hpull_B = frameTMP_B->pullHist("dataCTAUERR_Tot","pdfCTAUERR_Tot");
  //RooHist* hpull_B = frameTMP_B->pullHist("dataHist_Sig","pdfCTAUERR_Jpsi");
  hpull_B->SetMarkerSize(0.8);
  //RooPlot* pullFrame_B = ws->var("ctau3DErr")->frame(Title("Pull Distribution")) ;
  RooPlot* pullFrame_B = ws->var("ctau3DErr")->frame(Bins(nBins),Range(minRange-0.01, maxRange+0.01)) ;
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

  TLine *lB = new TLine(minRange-0.01,0, maxRange+0.01,0);
  lB->SetLineStyle(1);
  lB->Draw("same");
  pad_B_2->Update();
  cout << endl << "************** Finished SPLOT *****************" << endl << endl;

  c_B->Update();
  c_B->SaveAs(Form("figs/2DFit_%s/CtauErr_v2Bins_n/ctauErr_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_v2_%.2f.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP, v2));

  TFile *outFile = new TFile(Form("roots/2DFit_%s/CtauErr_v2Bins_n/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_v2_%.2f.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP, v2),"recreate");
  dataw_Bkg->Write();
  dataw_Sig->Write();
  pdfCTAUERR_Tot->Write();
  pdfCTAUERR_Bkg->Write();
  pdfCTAUERR_Jpsi->Write();
  outFile->Close();


}
