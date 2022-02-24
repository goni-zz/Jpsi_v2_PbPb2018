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

void CtauBkg_v2Bin_LowPt_x(
    double v2=0,
    double ptLow=3, double ptHigh=4.5,
    float yLow=1.6, float yHigh=2.4,
    int cLow=20, int cHigh=120,
    int PRw=1, bool fEffW = true, bool fAccW = true, bool isPtW = true, bool isTnP = true
    )
{


  TString DATE;
  if(ptLow==6.5&&ptHigh==50&&!(cLow==0&&cHigh==180)) DATE=Form("%i_%i",0,180);
  else DATE=Form("%i_%i",cLow/2,cHigh/2);
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir(Form("roots/2DFit_%s/CtauBkg_v2Bins_x",DATE.Data()),kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/CtauBkg_v2Bins_x",DATE.Data()),kTRUE);

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

  TFile* f1; TFile* fMass; TFile* fCErr; TFile* fCRes;
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);

  fMass = new TFile(Form("roots/2DFit_%s/Mass_v2Bins_x/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_v2_%.2f.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP, v2));
  fCErr = new TFile(Form("roots/2DFit_%s/CtauErr_v2Bins_x/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_v2_%.2f.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP, v2));
  fCRes = new TFile(Form("roots/2DFit_%s/CtauRes_v2Bins_x/CtauResResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_v2_%.2f.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP, v2));

  RooDataSet *datasetMass = (RooDataSet*)fMass->Get("datasetMass");
  RooAddPdf* pdfMASS_Tot = (RooAddPdf*)fMass->Get("pdfMASS_Tot");
  RooDataSet *dataw_Bkg = (RooDataSet*)fCErr->Get("dataw_Bkg");
  RooAddPdf* GaussModel_Tot = (RooAddPdf*)fCRes->Get("GaussModel_Tot");
  RooAddPdf* pdfCTAUERR_Bkg = (RooAddPdf*)fCErr->Get("pdfCTAUERR_Bkg");

  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*datasetMass);
  ws->import(*pdfMASS_Tot);
  ws->import(*dataw_Bkg);
  ws->import(*GaussModel_Tot);
  ws->import(*pdfCTAUERR_Bkg);
  double ctauErrMin, ctauErrMax;
  ws->data("dataw_Bkg")->getRange(*ws->var("ctau3DErr"), ctauErrMin, ctauErrMax);
  cout<<"####################################" <<endl;
  cout<<"pt: "<<ptLow<<"-"<<ptHigh<<", y: "<<yLow<<"-"<<yHigh<<", Cent: "<<cLow/2<<"-"<<cHigh/2<<"%"<<endl;

  cout <<"******** New Combined Dataset ***********" <<endl;
  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  ws->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);
  ws->var("ctau3DErr")->setRange(ctauErrMin, ctauErrMax);
  ws->var("ctau3DErr")->setRange("ctauErrRange",ctauErrMin, ctauErrMax);
  ws->var("ctau3D")->Print();
  ws->var("ctau3DErr")->Print();

  //***********************************************************************
  //**************************** Bkg CTAU FIT *****************************
  //***********************************************************************
  cout <<endl <<"************** Start BKG Ctau Fit *****************" <<endl <<endl;
  //make parameter 3 exp
  ws->factory("zeroMean[0.0]");
  //ws->factory("b_Bkg[0.1, 0., 1.]");//NP fraction for bkg
  //ws->factory("fDFSS[0.5, 1e-6, 1.]");
  //ws->factory("fDLIV[0.5, 1e-6, 1.]");
  //ws->factory("lambdaDDS_Bkg[0.2, 1e-6, 1.]");
  //ws->factory("lambdaDF_Bkg[0.3, 1e-6, 1.]");
  //ws->factory("lambdaDSS_Bkg[0.5, 1e-6, 1.]");
  //if(ptLow==3&&ptHigh==4.5){
  //ws->factory("b_Bkg[0.9480, 0., 1.]");//NP fraction for bkg
  //ws->factory("fDFSS[0.5909, 0., 1.]");
  //ws->factory("fDLIV[0.5640, 0., 1.]");
  //ws->factory("lambdaDDS_Bkg[0.2990, 1e-6, 1.]");
  //ws->factory("lambdaDF_Bkg1[0.0371, 1e-6, 1.]");
  //ws->factory("lambdaDF_Bkg2[1e-6, 1e-6, 1.]");
  //ws->factory("lambdaDSS_Bkg1[0.0591, 1e-6, 1.]");
  //ws->factory("lambdaDSS_Bkg2[0.3689, 1e-6, 1.]");
  //ws->factory("fDSS12[0.3, 0., 1.]");
  //ws->factory("fDF12[0.3, 0., 1.]");}
  //if(ptLow==3&&ptHigh==4.5){
  //ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  //ws->factory("fDFSS[0.6, 0., 1.]");
  //ws->factory("fDLIV[0.4, 0., 1]");
  //ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  //ws->factory("lambdaDF_Bkg1[0.0006, 1e-6, 1]");
  //ws->factory("lambdaDF_Bkg2[5e-3, 1e-6, 1]");
  //ws->factory("lambdaDSS_Bkg1[0.07, 1e-6, 1.]");
  //ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1.]");
  //ws->factory("fDSS12[0.5, 0., 1.]");
  //ws->factory("fDF12[0.5, 0., 1.]");}

//  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
//  ws->factory("fDFSS[0.6, 0., 1.]");
//  ws->factory("fDLIV[0.6, 0., 0.99]");
//  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
//  ws->factory("lambdaDF_Bkg1[0.006, 1e-6, 1.]");
//  ws->factory("lambdaDF_Bkg2[3e-5, 1e-6, 1.]");
//  ws->factory("lambdaDSS_Bkg1[0.1, 1e-6, 1.]");
//  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1.]");
//  ws->factory("fDSS12[0.5, 0., 1.]");
//  ws->factory("fDF12[0.5, 0., 1.]");
  if(v2==2.7&&ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.4, 0., 1]");
  ws->factory("lambdaDDS_Bkg[0.04, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.0327, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.1, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.3, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1.]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==-1.5&&ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.9480, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.89909, 0., 1.]");
  ws->factory("fDLIV[0.8640, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.2990, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.0371, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg2[1e-6, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg1[0.0591, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg2[0.3689, 1e-6, 1.]");
  ws->factory("fDSS12[0.3, 0., 1.]");
  ws->factory("fDF12[0.3, 0., 1.]");}
  else if(v2==0.&&ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.4, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.06, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg2[0.01, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg1[0.1, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.3, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==0.3&&ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.3, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.006, 1e-6, 0.1]");
  ws->factory("lambdaDF_Bkg2[5e-3, 1e-6, 0.1]");
  ws->factory("lambdaDSS_Bkg1[0.1, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==-0.3&&ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.3, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.06, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg2[0.01, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg1[0.5, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1.]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==-0.6&&ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.6, 0., 0.99]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.0006, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg2[3e-5, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg1[0.1, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1.]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==-0.9&&ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.6, 0., 0.99]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.0006, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg2[3e-5, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg1[0.1, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1.]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==0.9&&ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.6, 0., 0.99]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.0006, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg2[3e-5, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg1[0.1, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1.]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==-1.2&&ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.6, 0., 1]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.006, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg2[3e-5, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg1[0.1, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1.]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==-2.1&&ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.6, 0., 0.99]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.0006, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg2[3e-5, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg1[0.1, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1.]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==2.1&&ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.6, 0., 0.99]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.006, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg2[0.01, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg1[0.1, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1.]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==-0.3&&ptLow==4.5&&ptHigh==6.5){
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.3, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.016, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.1, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.2, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
//  else if(ptLow==4.5&&ptHigh==6.5){
//  ws->factory("b_Bkg[0.8, 0., 1.]");//NP fraction for bkg
//  ws->factory("fDFSS[0.6, 0., 1.]");
//  ws->factory("fDLIV[0.8, 0., 1.]");
//  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
//  ws->factory("lambdaDF_Bkg1[0.01, 1e-6, 1]");
//  ws->factory("lambdaDF_Bkg2[0.02, 1e-6, 1]");
//  ws->factory("lambdaDSS_Bkg1[0.2, 1e-6, 1]");
//  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1]");
//  ws->factory("fDSS12[0.5, 0., 1.]");
//  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==0.6&&ptLow==3&&ptHigh==6.5){
  cout << "++++++++++++++++++++++ HERE +++++++++++++++++++++++" << endl;
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.3, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.06, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.01, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.2, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==0.9&&ptLow==3&&ptHigh==6.5){
  cout << "++++++++++++++++++++++ HERE +++++++++++++++++++++++" << endl;
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.3, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.4, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.02, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.01, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.2, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==-0.6&&ptLow==4.5&&ptHigh==6.5){
  cout << "++++++++++++++++++++++ HERE +++++++++++++++++++++++" << endl;
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.3, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.06, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.01, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.2, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==0.9&&ptLow==4.5&&ptHigh==6.5){
  ws->factory("b_Bkg[0.8, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.8, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.04, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.11, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.2, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==1.8&&ptLow==4.5&&ptHigh==6.5){ //-1.8
  ws->factory("b_Bkg[0.8, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.8, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.2, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.01, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.02, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.3, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.05, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==2.1&&ptLow==4.5&&ptHigh==6.5){ //-1.8
  ws->factory("b_Bkg[0.8, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.8, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.2, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.01, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.01, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.3, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.05, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==-0.9&&ptLow==4.5&&ptHigh==6.5){
  ws->factory("b_Bkg[0.8, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.8, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.01, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.02, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.2, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==-1.5&&ptLow==4.5&&ptHigh==6.5){
  ws->factory("b_Bkg[0.8, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.8, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.01, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.02, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.2, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==-2.1&&ptLow==4.5&&ptHigh==6.5){ //-1.8
  ws->factory("b_Bkg[0.8, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.8, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.2, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.01, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.01, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.3, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.05, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==-2.4&&ptLow==4.5&&ptHigh==6.5){ //-1.8
  ws->factory("b_Bkg[0.8, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.7, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.03, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.03, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.2, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.05, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else if(v2==1.2&&ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.3, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.01, 1e-6, 1]");
  ws->factory("lambdaDF_Bkg2[0.01, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg1[0.2, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.4, 1e-6, 1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}
  else {
  ws->factory("b_Bkg[0.2, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.6, 0., 1.]");
  ws->factory("fDLIV[0.3, 0., 1.]");
  ws->factory("lambdaDDS_Bkg[0.1, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg1[0.006, 1e-6, 0.1]");
  ws->factory("lambdaDF_Bkg2[0.006, 1e-6, 0.1]");
  ws->factory("lambdaDSS_Bkg1[0.3, 1e-6, 1]");
  ws->factory("lambdaDSS_Bkg2[0.1, 1e-6, .1]");
  ws->factory("fDSS12[0.5, 0., 1.]");
  ws->factory("fDF12[0.5, 0., 1.]");}

//  ws->factory("b_Bkg[0.9480, 0., 1.]");//NP fraction for bkg
//  ws->factory("fDFSS[0.5909, 0., 1.]");
//  ws->factory("fDLIV[0.9640, 0., 1.]");
//  ws->factory("lambdaDDS_Bkg[0.2990, 1e-6, 1.]");
//  ws->factory("lambdaDF_Bkg1[0.0371, 1e-6, 1.]");
//  ws->factory("lambdaDF_Bkg2[1e-6, 1e-6, 1.]");
//  ws->factory("lambdaDSS_Bkg1[0.0591, 1e-6, 1.]");
//  ws->factory("lambdaDSS_Bkg2[0.3689, 1e-6, 1.]");
//  ws->factory("fDSS12[0.3, 0., 1.]");
//  ws->factory("fDF12[0.3, 0., 1.]");}
  //else {
  //ws->factory("b_Bkg[0.9480, 0., 1.]");//NP fraction for bkg
  //ws->factory("fDFSS[0.5909, 0., 1.]");
  //ws->factory("fDLIV[0.9640, 0., 1.]");
  //ws->factory("lambdaDDS_Bkg[0.2990, 1e-6, 1.]");
  //ws->factory("lambdaDF_Bkg1[0.0371, 1e-6, 1.]");
  //ws->factory("lambdaDF_Bkg2[1e-6, 1e-6, 1.]");
  //ws->factory("lambdaDSS_Bkg1[0.0591, 1e-6, 1.]");
  //ws->factory("lambdaDSS_Bkg2[0.3689, 1e-6, 1.]");
  //ws->factory("fDSS12[0.3, 0., 1.]");
  //ws->factory("fDF12[0.3, 0., 1.]");}
  //parameters fixed by Resolution model
  ws->var("ctau1_CtauRes")->setConstant(kTRUE); ws->var("s1_CtauRes")->setConstant(kTRUE);
  ws->var("ctau2_CtauRes")->setConstant(kTRUE);	ws->var("rS21_CtauRes")->setConstant(kTRUE);
  //ws->var("ctau3_CtauRes")->setConstant(kTRUE);	ws->var("rS32_CtauRes")->setConstant(kTRUE);
  //ws->var("f2_CtauRes")->setConstant(kTRUE);	
  ws->var("f_CtauRes")->setConstant(kTRUE);
  cout<<ws->var("N_Bkg")->getVal()<<"+/-"<<ws->var("N_Bkg")->getError()<<endl;
  cout<<ws->var("ctau1_CtauRes")->getVal()<<endl;
  cout<<ws->var("f_CtauRes")->getVal()<<"+/-"<<ws->var("f_CtauRes")->getError()<<endl;
  cout<<ws->var("s1_CtauRes")->getVal()<<endl;
  //make res model
  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes1", "ctau3D",
        "ctau1_CtauRes", //"ctau1_CtauRes",
        "s1_CtauRes",
        "zeroMean",
        "ctau3DErr"
        ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes2", "ctau3D",
        "ctau2_CtauRes", //"ctau2_CtauRes",
        "s2_CtauRes",
        "zeroMean",
        "ctau3DErr"
        ));
  //ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes3", "ctau3D",
  //      "ctau3_CtauRes", //"ctau3_CtauRes",
  //      "s3_CtauRes",
  //      "zeroMean",
  //      "ctau3DErr"
  //      ));
  //ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "ctauRes32", "ctauRes3", "ctauRes2", "f2_CtauRes"));
  ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "pdfCTAURES", "ctauRes1", "ctauRes2", "f_CtauRes"));
  //make 3 exp
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS1", "ctau3D", "lambdaDSS_Bkg1", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS2", "ctau3D", "lambdaDSS_Bkg2", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF1", "ctau3D", "lambdaDF_Bkg1", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF2", "ctau3D", "lambdaDF_Bkg2", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", "pdfCTAUDDS", "ctau3D", "lambdaDDS_Bkg", "pdfCTAURES"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUDSS", "fDSS12", "pdfCTAUDSS1", "pdfCTAUDSS2"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUDF", "fDF12", "pdfCTAUDF1", "pdfCTAUDF2"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU1", "fDFSS", "pdfCTAUDSS", "pdfCTAUDF"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_BkgNoPR", "fDLIV", "pdfCTAU1", "pdfCTAUDDS"));//NP
  ws->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_BkgPR","pdfCTAURES"));//PR
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_Bkg", "b_Bkg", "pdfCTAUCOND_BkgNoPR", "pdfCTAUCOND_BkgPR"));
  ws->factory(Form("RooExtendPdf::%s(%s,%s)", "pdfTot_Bkg", "pdfCTAUCOND_Bkg", "N_Bkg"));//N_Bkg is number of bkg from dataw_Bkg

  RooProdPdf pdfPR("pdfCTAU_BkgPR", "", *ws->pdf("pdfCTAUERR_Bkg"), Conditional( *ws->pdf("pdfCTAUCOND_BkgPR"), RooArgList(*ws->var("ctau3D")) ));
  ws->import(pdfPR);
  RooProdPdf pdfNoPR("pdfCTAU_BkgNoPR", "", *ws->pdf("pdfCTAUERR_Bkg"), Conditional( *ws->pdf("pdfCTAUCOND_BkgNoPR"), RooArgList(*ws->var("ctau3D")) ));
  ws->import(pdfNoPR);
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU_Bkg", "b_Bkg", "pdfCTAU_BkgNoPR", "pdfCTAU_BkgPR"));
  RooAbsPdf *pdfCTAU_Bkg_Tot = new RooAddPdf("pdfCTAU_Bkg_Tot", "pdfCTAU_Bkg_Tot", RooArgList(*ws->pdf("pdfCTAU_Bkg")), RooArgList(*ws->var("N_Bkg")));
  ws->import(*pdfCTAU_Bkg_Tot);
  //RooAbsPdf* ctauBkgModel = ctauBkgModel = new RooAddPdf("pdfTot_Bkg","pdfTot_Bkg",*ws->pdf("pdfCTAUCOND_Bkg"), *ws->var("N_Bkg"));
  //ws->import(*pdfCTAUCOND_Bkg);

  TH1D* hTot = (TH1D*)ws->data("dataw_Bkg")->createHistogram(("hTot"), *ws->var("ctau3D"),Binning(nCtauBins,ctauLow,ctauHigh));
  double ctauMin=hTot->GetBinLowEdge(hTot->FindFirstBinAbove(1,1));
  if (cLow==80&&cHigh==100) ctauMin=-1.0;
  else if (cLow==100&&cHigh==180) ctauMin=-0.6;
  double ctauMax=hTot->GetBinLowEdge(hTot->FindLastBinAbove(2,1))+hTot->GetBinWidth(hTot->FindLastBinAbove(2,1));
  if (v2==-1.2&&ptLow==3&&ptHigh==4.5) { ctauMin=-1.1; }
  else if (v2==2.1&&ptLow==3&&ptHigh==4.5) { ctauMin=-1; ctauMax=1;}
  else if(v2==-0.9&&ptLow==4.5&&ptHigh==6.5){ ctauMax=2.3; }
  else if(v2==-1.5&&ptLow==4.5&&ptHigh==6.5){ ctauMax=2.2; }
  else if(v2==-1.8&&ptLow==4.5&&ptHigh==6.5){ ctauMax=1.2; }
  else if(v2==-2.1&&ptLow==4.5&&ptHigh==6.5){ ctauMax=1.; }
  else if(v2==-2.4&&ptLow==4.5&&ptHigh==6.5){ ctauMin=-1.; ctauMax=1.; }
  else if(v2==0.9&&ptLow==4.5&&ptHigh==6.5){ ctauMax=2.; }
  else if(v2==1.8&&ptLow==4.5&&ptHigh==6.5){ ctauMax=2.5; }
  else if(v2==2.1&&ptLow==4.5&&ptHigh==6.5){ ctauMin=-1.; ctauMax=1.; }
  else if(v2==2.7&&ptLow==3&&ptHigh==4.5) {ctauMin=-1.; ctauMax=1.;}
  else if(v2==-2.1&&ptLow==3&&ptHigh==6.5) {ctauMin=-1.; ctauMax=2;}
  else if(v2==-2.7&&ptLow==3&&ptHigh==6.5) {ctauMin=-1.;}
//  if (v2==0.3&&ptLow==3&&ptHigh==4.5) { ctauMax=.8; ctauMin=-.5; }
//  if (v2==-2.4&&ptLow==3&&ptHigh==4.5) { ctauMax=1.0; ctauMin=-1.; }
//  if (v2==1.8&&ptLow==3&&ptHigh==4.5) { ctauMax=1.0; ctauMin=-1.; }
//  if (v2==1.2&&ptLow==4.5&&ptHigh==6.5) { ctauMax=2.0; }
//  if (v2==2.4&&ptLow==4.5&&ptHigh==6.5) { ctauMax=2.0; ctauMin=-1.;}


  TCanvas* c_E =  new TCanvas("canvas_E","My plots",1108,4,550,520);
  c_E->cd();
  TPad *pad_E_1 = new TPad("pad_E_1", "pad_E_1", 0, 0.16, 0.98, 1.0);
  pad_E_1->SetTicks(1,1);
  pad_E_1->Draw(); pad_E_1->cd();
  RooPlot* myPlot_E = ws->var("ctau3D")->frame(Bins(nCtauBins), Range(ctauLow, ctauHigh)); // bins
  myPlot_E->SetTitle("");

  ws->pdf("pdfCTAU_Bkg_Tot")->setNormRange("ctauWindow");

  //RooDataSet* dataToFit = (RooDataSet*)dataw_Bkg->reduce(Form("ctau3D>=%.f&&ctau3D<=%.f",ctauMin, ctauMax))->Clone("dataw_Bkg");
  RooDataSet* dataToFit = (RooDataSet*)dataw_Bkg->reduce(Form("ctau3D>=%.f&&ctau3D<=%.f",ctauMin, ctauMax))->Clone("dataw_Bkg");
  ws->import(*dataToFit, Rename("dataToFit"));

  pad_E_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot2_E = (RooPlot*)myPlot_E->Clone();

  double normDSTot = 1.0;
  if (ws->data("dataToFit")){normDSTot = ws->data("dataToFit")->sumEntries()/ws->data("dataToFit")->sumEntries();}
  double normBkg = 1.0;
  if (ws->data("dataw_Bkg")){normBkg = ws->data("dataToFit")->sumEntries()*normDSTot/ws->data("dataw_Bkg")->sumEntries();}

  bool isWeighted = ws->data("dataw_Bkg")->isWeighted();
  //RooFitResult* fitCtauBkg = ws->pdf("pdfTot_Bkg")->fitTo(*dataw_Bkg, Save(), Range("ctauRange"), Extended(kTRUE), NumCPU(nCPU), PrintLevel(-1));
  RooFitResult* fitCtauBkg = ws->pdf("pdfTot_Bkg")->fitTo(*dataToFit, Save(), Range("ctauRange"), Extended(kTRUE), NumCPU(nCPU), PrintLevel(-1), SumW2Error(isWeighted));
  ws->import(*fitCtauBkg, "fitCtauBkg");

  myPlot2_E->updateNormVars(RooArgSet(*ws->var("mass"), *ws->var("ctau3D"), *ws->var("ctau3DErr"))) ;

  ws->data("dataToFit")->plotOn(myPlot2_E, Name("data_ctauBkg"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlue), LineColor(kBlue), MarkerSize(0.7));
  ws->data("dataw_Bkg")->plotOn(myPlot2_E, Name("data_ctauBkg"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.7));
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(myPlot2_E, Name("ctauBkg_Tot"),  Normalization(normBkg, RooAbsReal::NumEvent), NumCPU(nCPU),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE),
      FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), Precision(1e-4));
  ws->data("dataToFit")->plotOn(myPlot2_E, Name("data_ctauBkg"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.7));
  if (ws->pdf("pdfCTAU_BkgPR")) {ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(myPlot2_E,Name("BKGPR"),Components(RooArgSet(*ws->pdf("pdfCTAU_BkgPR"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE), Normalization(normBkg, RooAbsReal::NumEvent), LineColor(kRed+2), Precision(1e-4), NumCPU(nCPU));}
  if (ws->pdf("pdfCTAU_BkgNoPR")) {ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(myPlot2_E,Name("BKGNoPR"),Components(RooArgSet(*ws->pdf("pdfCTAU_BkgNoPR"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE), Normalization(normBkg, RooAbsReal::NumEvent), LineColor(kOrange+10), Precision(1e-4), NumCPU(nCPU));}
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(myPlot2_E,Name("PDF"),  Normalization(normBkg, RooAbsReal::NumEvent), NumCPU(nCPU), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE), LineColor(kBlack), Precision(1e-4));

  myPlot2_E->GetYaxis()->SetRangeUser(10e-2, 10e7);
  TH1* h = ws->data("dataw_Bkg")->createHistogram("hist", *ws->var("ctau3D"), Binning(myPlot_E->GetNbinsX(),myPlot_E->GetXaxis()->GetXmin(),myPlot_E->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));
  Double_t Yup(0.),Ydown(0.);
  //Yup = YMax*TMath::Power((YMax/0.1), 0.5);
  Yup = YMax*TMath::Power((YMax/0.01), 0.5);
  Ydown = 0.01;
  myPlot2_E->GetYaxis()->SetRangeUser(Ydown,Yup);
  myPlot2_E->GetXaxis()->SetRangeUser(-4, 7);
  myPlot2_E->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  myPlot2_E->SetFillStyle(4000);
  myPlot2_E->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_E->GetXaxis()->SetLabelSize(0);
  myPlot2_E->GetXaxis()->SetTitleSize(0);

  TLine   *minline = new TLine(ctauMin, 0.0, ctauMin, Ydown*TMath::Power((Yup/Ydown),0.4));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  myPlot2_E->addObject(minline);
  TLine   *maxline = new TLine(ctauMax, 0.0, ctauMax, Ydown*TMath::Power((Yup/Ydown),0.4));
  maxline->SetLineStyle(2);
  maxline->SetLineColor(1);
  maxline->SetLineWidth(3);
  myPlot2_E->addObject(maxline);

  myPlot2_E->Draw();
  TLegend* leg_E = new TLegend(text_x+0.25,text_y+0.04,text_x+0.38,text_y-0.15); leg_E->SetTextSize(text_size);
  leg_E->SetTextFont(43);
  leg_E->SetBorderSize(0);
  leg_E->AddEntry(myPlot2_E->findObject("data_ctauBkg"),"Data_Bkg","pe");
  leg_E->AddEntry(myPlot2_E->findObject("ctauBkg_Tot"),"Total PDF","fl");
  //leg_E->AddEntry(myPlot2_E->findObject("test"),"?? PDF","l");
  leg_E->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);

  drawText(Form("N_{Bkg} = %.f #pm %.f", ws->var("N_Bkg")->getVal(), ws->var("N_Bkg")->getError() ),text_x+0.5,text_y,text_color,text_size);
  drawText(Form("b_{Bkg} = %.4f #pm %.4f", ws->var("b_Bkg")->getVal(), ws->var("b_Bkg")->getError() ),text_x+0.5,text_y-y_diff*1,text_color,text_size);
  drawText(Form("fDFSS = %.4f #pm %.4f", ws->var("fDFSS")->getVal(), ws->var("fDFSS")->getError() ),text_x+0.5,text_y-y_diff*2,text_color,text_size);
  drawText(Form("fDLIV = %.4f #pm %.4f", ws->var("fDLIV")->getVal(), ws->var("fDLIV")->getError() ),text_x+0.5,text_y-y_diff*3,text_color,text_size);
  drawText(Form("#lambdaDDS_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDDS_Bkg")->getVal(), ws->var("lambdaDDS_Bkg")->getError() ),text_x+0.5,text_y-y_diff*4,text_color,text_size);
  drawText(Form("#lambdaDF1_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDF_Bkg1")->getVal(), ws->var("lambdaDF_Bkg1")->getError() ),text_x+0.5,text_y-y_diff*5,text_color,text_size);
  drawText(Form("#lambdaDF2_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDF_Bkg2")->getVal(), ws->var("lambdaDF_Bkg2")->getError() ),text_x+0.5,text_y-y_diff*6,text_color,text_size);
  drawText(Form("#lambdaDSS1_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDSS_Bkg1")->getVal(), ws->var("lambdaDSS_Bkg1")->getError() ),text_x+0.5,text_y-y_diff*7,text_color,text_size);
  drawText(Form("#lambdaDSS2_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDSS_Bkg2")->getVal(), ws->var("lambdaDSS_Bkg2")->getError() ),text_x+0.5,text_y-y_diff*8,text_color,text_size);
  //pullDist(ws, pad_E_2, c_E, frameTMP_E, hpull_E, "data_ctauBkg", "ctauBkg_Tot", "ctau3D", nCtauBins, ctauLow, ctauHigh, "#font[12]{l}_{J/#psi} (mm)");
  
  TPad *pad_E_2 = new TPad("pad_E_2", "pad_E_2", 0, 0.006, 0.98, 0.227);
  c_E->cd();
  pad_E_2->Draw();
  pad_E_2->cd();
  pad_E_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_E_2->SetBottomMargin(0.67);
  pad_E_2->SetBottomMargin(0.4);
  pad_E_2->SetFillStyle(4000);
  pad_E_2->SetFrameFillStyle(4000);
  pad_E_2->SetTicks(1,1);

  RooPlot* frameTMP = (RooPlot*)myPlot2_E->Clone("TMP");
  RooHist* hpull_E = frameTMP->pullHist("data_ctauBkg","ctauBkg_Tot", true);
  hpull_E->SetMarkerSize(0.8);
  RooPlot* pullFrame_E = ws->var("ctau3D")->frame(Title("Pull Distribution"), Bins(nCtauBins), Range(ctauLow,ctauHigh)) ;
  pullFrame_E->addPlotable(hpull_E,"PX") ;
  pullFrame_E->SetTitle("");
  pullFrame_E->SetTitleSize(0);
  pullFrame_E->GetYaxis()->SetTitleOffset(0.3) ;
  pullFrame_E->GetYaxis()->SetTitle("Pull") ;
  pullFrame_E->GetYaxis()->SetTitleSize(0.15) ;
  pullFrame_E->GetYaxis()->SetLabelSize(0.15) ;
  pullFrame_E->GetYaxis()->SetRangeUser(-7,7);
  pullFrame_E->GetYaxis()->CenterTitle();
  pullFrame_E->GetYaxis()->SetTickSize(0.04);
  pullFrame_E->GetYaxis()->SetNdivisions(404);

  pullFrame_E->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  //pullFrame_E->GetXaxis()->SetRangeUser(-1, 7);
  pullFrame_E->GetXaxis()->SetTitleOffset(1.05) ;
  pullFrame_E->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame_E->GetXaxis()->SetLabelSize(0.15) ;
  pullFrame_E->GetXaxis()->SetTitleSize(0.15) ;
  pullFrame_E->GetXaxis()->CenterTitle();
  pullFrame_E->GetXaxis()->SetTickSize(0.03);
  pullFrame_E->Draw() ;

  TLine *lD = new TLine(ctauLow, 0, ctauHigh, 0);
  lD->SetLineStyle(1);
  lD->Draw("same");

  printChi2(ws, pad_E_2, frameTMP, fitCtauBkg, "ctau3D", "data_ctauBkg", "ctauBkg_Tot", nCtauBins, false);

  pad_E_2->Update();

  c_E->Update();
  c_E->SaveAs(Form("figs/2DFit_%s/CtauBkg_v2Bins_x/Bkg_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_v2_%.2f.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP, v2));
  RooArgSet* fitargs = new RooArgSet();
  fitargs->add(fitCtauBkg->floatParsFinal());
  RooDataSet *datasetCBkg = new RooDataSet("datasetCBkg","dataset with Ctau Background Fit result", *fitargs);
  RooWorkspace *wscbkg = new RooWorkspace("workspaceCBkg");
  wscbkg->import(*fitCtauBkg);
  //wscbkg->import(*pdfCTAUCOND_Bkg);

  //	ws->Print();

  TFile *outFile = new TFile(Form("roots/2DFit_%s/CtauBkg_v2Bins_x/CtauBkgResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_v2_%.2f.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP, v2),"recreate");
  fitCtauBkg->Write();
  pdfCTAU_Bkg_Tot->Write();
  //pdfCTAUCOND_Bkg->Write();
  //pdfTot_Bkg->Write();
  datasetCBkg->Write();
  wscbkg->Write();
  outFile->Close();
}
