#include <iostream>
#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
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
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit;

void Final2DFit_Corr_LowPt(
    float ptLow=3, float ptHigh=4.5,
    float yLow=1.6, float yHigh=2.4,
    int cLow=20, int cHigh=120,
    int PRw=1, bool fEffW = true, bool fAccW = true, bool isPtW = true, bool isTnP = true
    )
{

  TString DATE="10_60";
  //TString DATE="20_40";
  //TString DATE="0_180";
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir(Form("roots/2DFit_%s/Final",DATE.Data()),kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/Final", DATE.Data()),kTRUE);

  TString fname;
  if (PRw==1) fname="PR";
  else if (PRw==2) fname="NP";

  cout<<"pt: "<<ptLow<<" - "<<ptHigh<<", y: "<<yLow<<" - "<<yHigh<<", Cent: "<<cLow/2<<" - "<<cHigh/2<<"%"<<endl;

  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  TFile* f1; TFile* fMass; TFile* fCErr; TFile* fCRes; TFile* fCBkg; TFile* fCTrue;
  TString kineCut; TString OS;
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);

  //f1 = new TFile(Form("../../skimmedFiles/OniaRooDataSet_isMC0_JPsi_%sw_Effw%d_Accw%d_PtW%d_TnP%d_20210426.root",fname.Data(),fEffW,fAccW,isPtW,isTnP));
  f1 = new TFile("../../skimmedFiles/OniaRooDataSet_isMC0_JPsi_PRw_Effw1_Accw1_PtW1_TnP1_210915.root");
  fMass = new TFile(Form("roots/2DFit_%s/Mass/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  fCErr = new TFile(Form("roots/2DFit_%s/CtauErr/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  fCRes = new TFile(Form("roots/2DFit_%s/CtauRes/CtauResResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  fCBkg = new TFile(Form("roots/2DFit_%s/CtauBkg/CtauBkgResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  if(DATE=="0_180") fCTrue = new TFile(Form("roots/2DFit_%s/CtauTrue/CtauTrueResult_Inclusive_pt6.5-50.0_y0.0-2.4_muPt0.0_centrality0-180.root",DATE.Data()));
  else fCTrue = new TFile(Form("roots/2DFit_%s/CtauTrue/CtauTrueResult_Inclusive_%s.root",DATE.Data(),kineLabel.Data()));

  OS="recoQQsign==0 &&";
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut
  kineCut  = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.6 && cBin>%d && cBin<%d",ptLow, ptHigh, yLow, yHigh, cLow, cHigh);
  kineCut = OS+accCut+kineCut;

  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooDataSet *datasetMass = (RooDataSet*)fMass->Get("datasetMass");
  RooAddPdf* pdfMASS_Tot = (RooAddPdf*)fMass->Get("pdfMASS_Tot");
  RooDataSet *dataw_Bkg = (RooDataSet*)fCErr->Get("dataw_Bkg");
  //RooDataSet *dataw_Sig = (RooDataSet*)fCErr->Get("dataw_Sig");
  RooHistPdf* pdfCTAUERR_Tot = (RooHistPdf*)fCErr->Get("pdfCTAUERR_Tot");
  RooHistPdf* pdfCTAUERR_Jpsi = (RooHistPdf*)fCErr->Get("pdfCTAUERR_Jpsi");
  RooHistPdf* pdfCTAUERR_Bkg = (RooHistPdf*)fCErr->Get("pdfCTAUERR_Bkg");
  RooAddPdf* GaussModel_Tot = (RooAddPdf*)fCRes->Get("GaussModel_Tot");
  RooAddPdf* TrueModel_Tot = (RooAddPdf*)fCTrue->Get("TrueModel_Tot");
  RooAddPdf* pdfCTAU_Bkg_Tot = (RooAddPdf*)fCBkg->Get("pdfCTAU_Bkg_Tot");

  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataw_Bkg);
  double ctauErrMin;
  double ctauErrMax;
  
  //if (ptLow==7.5 && cLow==40){
  //ctauErrMin = 0.015;//0.0143889
  //ctauErrMax = 0.09;}
  //else if (ptLow==3 && ptHigh==6.5 && cLow==40){
  //ctauErrMin = 0.02;//0.0143889
  //ctauErrMax = 0.16;}
  //else if (ptLow==6 && ptHigh==7){
  //ctauErrMin = 0.021;//0.0143889
  //ctauErrMax = 0.09;}//0.124872
  //else if (ptLow==9 && cLow==40){
  //ctauErrMin = 0.015;//0.0143889
  //ctauErrMax = 0.08;}//0.124872
  //else if (ptLow==12 && cLow==40){
  //ctauErrMin = 0.01;//0.0143889
  //ctauErrMax = 0.05;}//0.124872
  //else if (ptLow==6.5 && cLow==40 && cHigh==60){
  //  ctauErrMin = ws->var("ctau3DErr")->getMin();
  //  ctauErrMax = 0.112428;}//0.124872
  //else {ctauErrMin = ws->var("ctau3DErr")->getMin();  ctauErrMax = ws->var("ctau3DErr")->getMax();}
  ctauErrMin = ws->var("ctau3DErr")->getMin();  ctauErrMax = ws->var("ctau3DErr")->getMax();
  //if (ptLow==3 && ptHigh==4.5 && cLow==20 && cHigh==120){
  //ctauErrMin = 0.03;//0.0143889
  //ctauErrMax = 0.16;}
  ws->import(*dataset); //total
  ws->import(*datasetMass);
  ws->import(*pdfMASS_Tot);
  //ws->import(*dataw_Sig);
  //ws->import(*GaussModel_Tot);
  ws->import(*TrueModel_Tot);
  ws->import(*pdfCTAUERR_Tot);
  ws->import(*pdfCTAUERR_Jpsi);
  //ws->import(*pdfCTAUERR_Bkg);
  ws->import(*pdfCTAU_Bkg_Tot);
  
  cout<<"CtauErr Min: "<<ctauErrMin<<", Max: "<<ctauErrMax<<endl;
  
  cout << "####################################" << endl;
  RooArgSet *argSet = new RooArgSet( *(ws->var("ctau3D")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("weight")), *(ws->var("ctau3DRes")), *(ws->var("ctau3DErr")) );
  argSet->add(*(ws->var("pt1")) ); argSet->add(*(ws->var("pt2")) ); argSet->add(*(ws->var("eta1")) );  argSet->add(*(ws->var("eta2")) ); argSet->add(*(ws->var("recoQQsign")) ); argSet->add(*(ws->var("cBin")) );
  RooDataSet *datasetW = new RooDataSet("datasetW","A sample",   *argSet, Import(*dataset),WeightVar(*ws->var("weight")));
  RooDataSet *datasetWo = new RooDataSet("datasetWo","A sample", *argSet, Import(*dataset));
  ws->import(*datasetW);
  ws->import(*datasetWo);
  
  RooDataSet *dsTot = (RooDataSet*)datasetW->reduce(*argSet, kineCut.Data() );
  //RooDataSet *dsTot = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data());
  //RooDataSet *dsToFit = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), Form("ctau3DErr>=%.6f&&ctau3DErr<=%.6f&&%s", ctauErrMin, ctauErrMax, kineCut.Data()) );
  //RooDataSet* dsToFit = (RooDataSet*)dsTot->reduce(Form("ctau3DErr>=%.6f&&ctau3DErr<=%.6f",0.009, 0.3))->Clone("dsTot");
  //RooDataSet* dsToFit = (RooDataSet*)dsTot->reduce(Form("ctau3DErr>=0.0126009&&ctau3DErr<=0.09"))->Clone("dsTot");
  dsTot->SetName("dsTot");
  ws->import(*dsTot);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  ws->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);
  //ws->var("ctau3DErr")->setRange(0.0126009, 0.09);
  ws->var("ctau3DErr")->setRange(ctauErrMin, ctauErrMax);
  //ws->var("ctau3DErr")->setRange("ctauRange", 0.0126009, 0.09);
  ws->var("ctau3DErr")->setRange("ctauRange", ctauErrMin, ctauErrMax);
  ws->var("ctau3DRes")->setRange(-10, 10);
  ws->var("ctau3DRes")->setRange("ctauRange", -10, 10);
  ws->var("mass")->Print();
  ws->var("ctau3D")->Print();
  ws->var("ctau3DErr")->Print();
  ws->var("ctau3DRes")->Print();
  //***********************************************************************
  //*************************** DRAW CTAU FIT *****************************
  //***********************************************************************
  //ws->pdf("pdfMASS_Tot")->getParameters(RooArgSet(*ws->var("mass")))->setAttribAll("Constant", kTRUE);
  ws->pdf("pdfMASS_Tot")->getParameters(
      RooArgSet(*ws->var("mass"), *ws->pdf("pdfMASS_Jpsi"), *ws->pdf("pdfMASS_bkg")
        ))->setAttribAll("Constant", kTRUE);
  //ws->pdf("GaussModelCOND_ctauRes")->getParameters(
  //    RooArgSet(*ws->var("ctau1_CtauRes"),*ws->var("ctau2_CtauRes"),*ws->var("ctau3_CtauRes"),
  //      *ws->var("s1_CtauRes"), *ws->var("rS21_CtauRes"), *ws->var("rS32_CtauRes"),
  //      *ws->var("f_CtauRes"), *ws->var("f2_CtauRes")
  //      ))->setAttribAll("Constant", kTRUE);
  ws->pdf("pdfCTAU_Bkg_Tot")->getParameters(
      //RooArgSet(*ws->var("ctau3D"), *ws->pdf("pdfCTAU_BkgPR"), *ws->pdf("pdfCTAU_BkgNoPR"),  *ws->pdf("pdfCTAUCOND_Bkg"), *ws->pdf("pdfCTAUCOND_BkgPR"), *ws->pdf("pdfCTAUCOND_BkgNoPR")
      RooArgSet(*ws->var("ctau3D"), *ws->pdf("pdfCTAU_BkgPR"), *ws->pdf("pdfCTAU_BkgNoPR"), *ws->pdf("pdfCTAUCOND_BkgPR"), *ws->pdf("pdfCTAUCOND_BkgNoPR")
        ))->setAttribAll("Constant", kTRUE);
  //ws->var("lambdaDSS1")->setConstant(kTRUE);//make it as a initial value..
  double lambda = ws->var("lambdaDSS")->getVal();
  double lambda1 = ws->var("lambdaDSS2")->getVal();
  double lambda2 = ws->var("lambdaDSS3")->getVal();
  double fdss = ws->var("fDSS")->getVal();
  double fdss1 = ws->var("fDSS1")->getVal();
  //double lambda2 = ws->var("lambdaDSS2")->getVal();
  //ws->var("b_Bkg")->setConstant(kTRUE);//
  //make jpsi pdf
  
  /*
  ws->factory(Form("lambdaDSS_test1[%.4f, %.4f, %.4f]", lambda, 1e-8, lambda*2));
  ws->factory(Form("lambdaDSS_test2[%.4f, %.4f, %.4f]", lambda1, 1e-8, lambda1*2));
  ws->factory(Form("lambdaDSS_test3[%.4f, %.4f, %.4f]", lambda2, 1e-8, lambda2*2));
  ws->factory(Form("fDSS1_test[%.4f, %.4f, %.4f]", fdss, 1e-8, 1.));
  ws->factory(Form("fDSS2_test[%.4f, %.4f, %.4f]", fdss1, 1e-8, 1.));
  */
  
  //ws->factory(Form("lambdaDSS_test1[%.4f, %.4f, %.4f]", lambda, lambda, lambda));
  //ws->factory(Form("lambdaDSS_test2[%.4f, %.4f, %.4f]", lambda1, lambda1, lambda1));
  //ws->factory(Form("lambdaDSS_test3[%.4f, %.4f, %.4f]", lambda2, lambda2, lambda2));
  //ws->factory(Form("fDSS1_test[%.4f, %.4f, %.4f]", fdss, fdss, fdss));
  //ws->factory(Form("fDSS2_test[%.4f, %.4f, %.4f]", fdss1, fdss1, fdss1));
  
  ws->factory(Form("lambdaDSS_test1[%.4f]", lambda));
  ws->factory(Form("lambdaDSS_test2[%.4f]", lambda1));
  ws->factory(Form("lambdaDSS_test3[%.4f]", lambda2));
  ws->factory(Form("fDSS1_test[%.4f]", fdss));
  ws->factory(Form("fDSS2_test[%.4f]", fdss1));

  //NoPR{
  //1exp
  //ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUCOND_JpsiNoPR", "ctau3D", "lambdaDSS", "pdfCTAURES")); //NP
  //3exp
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUE_test1", "ctau3D", "lambdaDSS_test1", "pdfCTAURES")); //NP
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUE_test2", "ctau3D", "lambdaDSS_test2", "pdfCTAURES")); //NP
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUE_test3", "ctau3D", "lambdaDSS_test3", "pdfCTAURES")); //NP
  ws->factory(Form("AddModel::%s({%s , %s}, %s)", "pdfCTAUTRUE_test12", "pdfCTAUTRUE_test1", "pdfCTAUTRUE_test2", "fDSS1_test"));
  ws->factory(Form("AddModel::%s({%s , %s}, %s)", "pdfCTAUCOND_JpsiNoPR", "pdfCTAUTRUE_test12", "pdfCTAUTRUE_test3", "fDSS2_test"));//}
  //PR
  ws->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_JpsiPR", "pdfCTAURES"));
  //3-4.5
  //if(ptLow==3&&ptHigh==4.5){ws->factory("b_Jpsi[0.09, 1e-8, 1]");}//NP fraction for Sig
  //else if(ptLow==4.5&&ptHigh==6.5){ws->factory("b_Jpsi[0.2, 1e-6, 1.]");}//NP fraction for Sig
  //else if(ptLow==6.5&&ptHigh==9){ws->factory("b_Jpsi[0.19, 0.19, .4]");}//NP fraction for Sig
  //else if(ptLow==6.5&&ptHigh==8.5&&yHigh==0.6){ws->factory("b_Jpsi[0.27, 0.21, .29]");}//NP fraction for Sig
  //else if(ptLow==6.5&&ptHigh==7.5)ws->factory("b_Jpsi[0.2487, 1e-8, .2988]");//NP fraction for Sig
  //else if(ptLow==7.5&&ptHigh==9){ws->factory("b_Jpsi[0.27, 0.21, .5]");}//NP fraction for Sig
  //else if(ptLow==9&&ptHigh==10){ws->factory("b_Jpsi[0.27, 0.0, 1.]");}//NP fraction for Sig
  //else if(ptLow==9&&ptHigh==12){ws->factory("b_Jpsi[0.29, 0.27, 0.43]");}//NP fraction for Sig
  //else if(ptLow==10&&ptHigh==12){ws->factory("b_Jpsi[0.3368, 0.25, 0.34]");}//NP fraction for Sig
  //if(ptLow==12&&ptHigh==15){ws->factory("b_Jpsi[0.4057, 0.38, 0.4725]");}//NP fraction for Sig
  //else if(ptLow==15&&ptHigh==50){ws->factory("b_Jpsi[0.48, 0.41, 0.53]");}//NP fraction for Sig
  //else if(ptLow==12&&ptHigh==50){ws->factory("b_Jpsi[0.49, 0.4, 0.6]");}//NP fraction for Sig
  //else ws->factory("b_Jpsi[0.1, 1e-8, 1]");//NP fraction for Sig
  if(ptLow==4&&ptHigh==5)ws->factory("b_Jpsi[0.20, 1e-3, 1.]");//NP fraction for Sig
  else if(ptLow==5&&ptHigh==6)ws->factory("b_Jpsi[0.18, 1e-3, 1.]");//NP fraction for Sig
  else ws->factory("b_Jpsi[0.31, 1e-3, 1.]");//NP fraction for Sig
  //ws->factory("b_Jpsi[0.05, 1e-8, 0.18]");//NP fraction for Sig

  //RooProdPdf pdfbkgPR("pdfCTAU_BkgPR", "", *ws->pdf("pdfCTAUERR_Bkg"),
  //    Conditional( *ws->pdf("pdfCTAUCOND_BkgPR"), RooArgList(*ws->var("ctau3D"))));
  //ws->import(pdfbkgPR);
  //RooProdPdf pdfbkgNoPR("pdfCTAU_BkgNoPR", "", *ws->pdf("pdfCTAUERR_Bkg"),
  //    Conditional( *ws->pdf("pdfCTAUCOND_BkgNoPR"), RooArgList(*ws->var("ctau3D"))));
  //ws->import(pdfbkgNoPR);

  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgPR",
        "pdfCTAU_BkgPR",
        "pdfMASS_bkg"
        ));
  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgNoPR",
        "pdfCTAU_BkgNoPR",
        "pdfMASS_bkg"
        ));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUMASS_Bkg",
        "b_Bkg",
        "pdfCTAUMASS_BkgNoPR",
        "pdfCTAUMASS_BkgPR"
        ));

  RooProdPdf pdfJpsiPR("pdfCTAU_JpsiPR","",*ws->pdf("pdfCTAUERR_Jpsi"),
      Conditional(*ws->pdf("pdfCTAUCOND_JpsiPR"),RooArgList(*ws->var("ctau3D"))));
  ws->import(pdfJpsiPR);
  RooProdPdf pdfJpsiNoPR("pdfCTAU_JpsiNoPR","",*ws->pdf("pdfCTAUERR_Jpsi"),
      Conditional(*ws->pdf("pdfCTAUCOND_JpsiNoPR"),RooArgList(*ws->var("ctau3D"))));
  ws->import(pdfJpsiNoPR);

  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_JpsiPR",
        "pdfCTAU_JpsiPR",
        "pdfMASS_Jpsi"
        ));
  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_JpsiNoPR",
        "pdfCTAU_JpsiNoPR",
        "pdfMASS_Jpsi"
        ));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUMASS_Jpsi",
        "b_Jpsi",
        "pdfCTAUMASS_JpsiNoPR",
        "pdfCTAUMASS_JpsiPR"
        ));
  RooAbsPdf *themodel =NULL;

  //double njpsi = ws->var("N_Jpsi")->getVal();
  //ws->factory(Form("N_Jpsi[%.3f, %.3f, %.3f]",njpsi, njpsi*0.9, njpsi*1.1));

  themodel = new RooAddPdf("pdfCTAUMASS_Tot", "pdfCTAUMASS_Tot",
      RooArgList(*ws->pdf("pdfCTAUMASS_Bkg"), *ws->pdf("pdfCTAUMASS_Jpsi")),
      RooArgList(*ws->var("N_Bkg"), *ws->var("N_Jpsi")) );
  ws->import(*themodel);
  //ws->pdf("pdfMASS_Tot")->getParameters(RooArgSet(*ws->var("mass")))->setAttribAll("Constant", kTRUE);
  std::vector< std::string > objs = {"Bkg", "Jpsi"};
  RooArgSet pdfList = RooArgSet("ConstraionPdfList");
  for (auto obj : objs) {
    if (ws->var(Form("N_%s", obj.c_str())))  {
      ws->factory(Form("Gaussian::%s_Gauss(%s,%s_Mean[%f],%s_Sigma[%f])",
            Form("N_%s", obj.c_str()), Form("N_%s", obj.c_str()),
            Form("N_%s", obj.c_str()), ws->var(Form("N_%s", obj.c_str()))->getValV(),
            Form("N_%s", obj.c_str()), ws->var(Form("N_%s", obj.c_str()))->getError()));

      pdfList.add(*ws->pdf(Form("N_%s_Gauss", obj.c_str())), kFALSE);
      std::cout << "[INFO] Constraining N_" << obj << " with Mean : " << ws->var(Form("N_%s_Mean", obj.c_str()))->getVal()
        << " and Sigma: " << ws->var(Form("N_%s_Sigma", obj.c_str()))->getVal() << std::endl;
    }
  }
  ws->defineSet("ConstrainPdfList", pdfList);
  ws->pdf("pdfCTAURES")->getParameters(RooArgSet(*ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("ctau3DErr")))->setAttribAll("Constant", kTRUE);
  cout<<ws->pdf("pdfCTAURES")->getVal()<<endl;
  ws->pdf("pdfCTAU_JpsiPR")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr")))->setAttribAll("Constant", kTRUE);
  ws->pdf("pdfCTAU_BkgPR")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr")))->setAttribAll("Constant", kTRUE);
  ws->pdf("pdfCTAU_BkgNoPR")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr")))->setAttribAll("Constant", kTRUE);

  RooArgSet* params = (RooArgSet*) ws->pdf("pdfCTAUMASS_Tot")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("mass"), *ws->var("ctau3DErr")));
  ws->saveSnapshot(("pdfCTAUMASS_Tot_parIni"),*params,kTRUE);
  delete params;

  RooArgSet *newpars = (RooArgSet*) ws->pdf("pdfCTAUMASS_Tot")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr"), *ws->var("ctau3DRes"), *ws->var("mass")));

  //TCanvas* c_G =  new TCanvas("canvas_G","My plots",880,800);
  TCanvas* c_G =  new TCanvas("canvas_G","My plots",1108,4,550,520);
  c_G->cd();
  TPad *pad_G_1 = new TPad("pad_G_1", "pad_G_1", 0, 0.16, 0.98, 1.0);
  pad_G_1->SetTicks(1,1);
  pad_G_1->Draw(); pad_G_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot_G = ws->var("ctau3D")->frame(nCtauBins*2); // bins
  //RooPlot* myPlot_G = ws->var("ctau3D")->frame(nCtauBins); // bins
  //RooPlot* myPlot_H = ws->var("mass")->frame(Bins(nMassBin), Range(2.6,3.5)); // bins
  myPlot_G->SetTitle("");

  c_G->cd();
  c_G->SetLogy();
  pad_G_1->cd();

  RooDataSet* dsToFit = (RooDataSet*)dsTot->reduce(Form("ctau3DErr>=%.6f&&ctau3DErr<=%.6f",ctauErrMin, ctauErrMax))->Clone("dsTot");
  dsToFit->SetName("dsToFit");
  ws->import(*dsToFit);
  //double normDSTot = ws->data("dsToFit")->sumEntries()/ws->data("dsTot")->sumEntries();
  //cout<<normDSTot<<": "<<ws->data("dsToFit")->sumEntries()<<"/"<<ws->data("dsTot")->sumEntries()<<endl;
  double normDSTot = (ws->var("N_Jpsi_Mean")->getVal()+ws->var("N_Bkg_Mean")->getVal())/ws->data("dsToFit")->sumEntries();
  cout<<normDSTot<<": "<<ws->var("N_Jpsi_Mean")->getVal()+ws->var("N_Bkg_Mean")->getVal()<<"/"<<ws->data("dsToFit")->sumEntries()<<endl;
  //double normBkg = ws->data("dsToFit")->sumEntries()*normDSTot/ws->data("dataw_Bkg")->sumEntries();
  //double normJpsi =ws->data("dsToFit")->sumEntries()*normDSTot/ws->data("dataw_Sig")->sumEntries();

  cout<<"##############START TOTAL CTAU FIT############"<<endl;
  RooFitResult* fitResult = ws->pdf("pdfCTAUMASS_Tot")->fitTo(*dsToFit, Extended(kTRUE), ExternalConstraints(*ws->set("ConstrainPdfList")),
      NumCPU(nCPU), SumW2Error(true), PrintLevel(3), Save());
  ws->import(*fitResult, "fitResult_pdfCTAUMASS_Tot");

  //DRAW
  RooPlot* myPlot2_G = (RooPlot*)myPlot_G->Clone();
  ws->data("dsToFit")->plotOn(myPlot_G,Name("dataHist_ctau"));
  myPlot2_G->updateNormVars(RooArgSet(*ws->var("mass"), *ws->var("ctau3D"), *ws->var("ctau3DErr")));
  //ws->data("dsToFit")->plotOn(myPlot2_G,Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G,Name("Ctau_Tot"),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE),
      Normalization(normDSTot, RooAbsReal::NumEvent),
      FillStyle(1001), FillColor(kViolet+6), VLines(), DrawOption("LF"), NumCPU(nCPU), LineColor(kBlack)
      );
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G,Name("Ctau_Bkg"),Components(RooArgSet( *ws->pdf("pdfCTAUMASS_Bkg") )),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE),
      Normalization(normDSTot, RooAbsReal::NumEvent),
      FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("F"), NumCPU(nCPU)
      );
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G,Name("Ctau_JpsiPR"),Components(RooArgSet( *ws->pdf("pdfCTAUMASS_JpsiPR") )),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE),
      Normalization(normDSTot, RooAbsReal::NumEvent),
      LineColor(kRed+3), NumCPU(nCPU)
      );
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G,Name("Ctau_JpsiNoPR"),Components(RooArgSet( *ws->pdf("pdfCTAUMASS_JpsiNoPR"))),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE),
      Normalization(normDSTot, RooAbsReal::NumEvent),
      LineColor(kGreen+3), NumCPU(nCPU)
      );
  ws->data("dsToFit")->plotOn(myPlot2_G,Name("data_Ctau"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack));
  //ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G,Name("Mass_TotLINE"),
  //                                     ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE),
  //                                     Normalization(normDSTot, RooAbsReal::NumEvent),
  //                                     LineColor(kBlack), NumCPU(32)
  //                                     );
  ws->saveSnapshot("pdfCTAUMASS_Tot_parFit",*newpars,kTRUE);
  myPlot2_G->GetYaxis()->SetRangeUser(10e-2, 10e8);
  TH1* hCtau = ws->data("dsToFit")->createHistogram("histCtau", *ws->var("ctau3D"), Binning(myPlot_G->GetNbinsX(),myPlot_G->GetXaxis()->GetXmin(),myPlot_G->GetXaxis()->GetXmax()));
  Double_t YMaxCtau = hCtau->GetBinContent(hCtau->GetMaximumBin());
  Double_t YMinCtau = 1e99;
  for (int i=1; i<=hCtau->GetNbinsX(); i++) if (hCtau->GetBinContent(i)>0) YMinCtau = min(YMinCtau, hCtau->GetBinContent(i));
  Double_t YupCtau(0.),YdownCtau(0.);
  YupCtau = YMaxCtau*TMath::Power((YMaxCtau/0.1), 0.5);
  YdownCtau = 0.1;
  myPlot2_G->GetYaxis()->SetRangeUser(YdownCtau,YupCtau);
  myPlot2_G->GetXaxis()->SetRangeUser(-4, 6);
  myPlot2_G->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  myPlot2_G->SetFillStyle(4000);
  myPlot2_G->GetXaxis()->SetLabelSize(0);
  myPlot2_G->GetXaxis()->SetTitleSize(0);
  myPlot2_G->Draw();
  TLegend* leg_G = new TLegend(text_x+0.25,text_y+0.05,text_x+0.38,text_y-0.11); leg_G->SetTextSize(text_size);
  leg_G->SetTextFont(43);
  leg_G->SetBorderSize(0);
  leg_G->AddEntry(myPlot2_G->findObject("data_Ctau"),"Data","pe");
  leg_G->AddEntry(myPlot2_G->findObject("Ctau_Tot"),"Total fit","fl");
  leg_G->AddEntry(myPlot2_G->findObject("Ctau_Bkg"),"Background","fl");
  leg_G->AddEntry(myPlot2_G->findObject("Ctau_JpsiPR"),"J/#psi Prompt","l");
  leg_G->AddEntry(myPlot2_G->findObject("Ctau_JpsiNoPR"),"J/#psi Non-Prompt","l");
  leg_G->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);
  drawText(Form("N_{J/#psi} = %.f #pm %.f",ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError()),text_x+0.5,text_y+0.05-y_diff,text_color,text_size);
  drawText(Form("N_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(), ws->var("N_Bkg")->getError() ),text_x+0.5,text_y+0.05-y_diff*2,text_color,text_size);
  drawText(Form("b_{J/#psi} = %.5f #pm %.5f",ws->var("b_Jpsi")->getVal(),ws->var("b_Jpsi")->getError()),text_x+0.5,text_y+0.05-y_diff*3,text_color,text_size);

  TPad *pad_G_2 = new TPad("pad_G_2", "pad_G_2", 0, 0.006, 0.98, 0.227);
  RooPlot* frameTMP_G = (RooPlot*)myPlot2_G->Clone("TMP_G");
  RooHist* hpull_G;
  pullDist(ws, pad_G_2, c_G, frameTMP_G, hpull_G, "data_Ctau", "Ctau_Tot", "ctau3D", nCtauBins, -4, 6.0, "#font[12]{l}_{J/#psi} (mm)");
  printChi2(ws, pad_G_2, frameTMP_G, fitResult, "ctau3D", "data_Ctau", "Ctau_Tot", nCtauBins);
  pad_G_2->Update();

  //TCanvas* c_H =  new TCanvas("canvas_H","My plots",880,800);
  TCanvas* c_H =  new TCanvas("canvas_H","My plots",1108,565,550,520);
  c_H->cd();
  TPad *pad_H_1 = new TPad("pad_H_1", "pad_H_1", 0, 0.16, 0.98, 1.0);
  pad_H_1->SetTicks(1,1);
  pad_H_1->Draw(); pad_H_1->cd();
  RooPlot* myPlot_H = ws->var("mass")->frame(Bins(nMassBin)); // bins
  myPlot_H->SetTitle("");
  c_H->cd();
  c_H->SetLogy();
  pad_H_1->cd();
  gPad->SetLogy();
  //RooAbsReal* fsigregion_model = pdfCTAUMASS_JpsiPR.createIntegral(*ws->var("ctau3D"),NormSet(*ws->var("ctau3D")),Range("ctauRange")); 
  //The "NormSet(x)" normalizes it to the total number of events to give you the fraction n_signal_region_events/n_total_events
  //RooArgSet* cloneSet = (RooArgSet*)RooArgSet(*ws->pdf("pdfCTAUMASS_JpsiPR"),"pdfCTAUMASS_JpsiPR").snapshot(kTRUE);
  //auto clone_mass_pdf = (RooAbsPdf*)cloneSet->find("pdfCTAUMASS_JpsiPR");

  //RooAbsReal* fsigregion_bkg = pdfCTAUMASS_JpsiNoPR.createIntegral(*ws->var("ctau3D"),NormSet(*ws->var("ctau3D")),Range("ctauRange"));
  //Double_t nsigevents = fsigregion_model->getVal()*(sig_yield.getVal()+bkg_yield.getVal())-fsigregion_bkg->getVal()*bkg_yield.getVal();
  //Double_t fsig = nsigevents/(fsigregion_model->getVal()*(sig_yield.getVal()+bkg_yield.getVal()));

  //themodel->RooAddPdf::pdfList().find("pdfCTAUMASS_JpsiPR");
  //RooAbsReal* PR_Frac = *ws->pdf("pdfCTAUMASS_JpsiPR")->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range(-4, 6.5));
  //RooAbsReal* NP_Frac = *ws->pdf("pdfCTAUMASS_JpsiNoPR")->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range(-4, 6.5));

  RooPlot* myPlot2_H = (RooPlot*)myPlot_H->Clone();
  ws->data("dsToFit")->plotOn(myPlot_H,Name("dataHist_mass"));
  myPlot2_H->updateNormVars(RooArgSet(*ws->var("mass"), *ws->var("ctau3D"), *ws->var("ctau3DErr")));

  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_H,Name("Mass_Bkg"),Components(RooArgSet( *ws->pdf("pdfCTAUMASS_Bkg") )),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE),
      Normalization(normDSTot, RooAbsReal::NumEvent),
      FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LF"), NumCPU(nCPU), LineStyle(kDashed)
      //FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), NumCPU(nCPU), LineStyle(kDashed)
      );
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_H,Name("Mass_JpsiPR"),Components(RooArgSet( *ws->pdf("pdfCTAUMASS_JpsiPR"), *ws->pdf("pdfCTAUMASS_Bkg") )),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE),
      Normalization(normDSTot, RooAbsReal::NumEvent),
      LineColor(kRed+3), LineStyle(1), Precision(1e-4), NumCPU(nCPU)
      );
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_H,Name("Mass_JpsiNoPR"),Components(RooArgSet( *ws->pdf("pdfCTAUMASS_JpsiNoPR"), *ws->pdf("pdfCTAUMASS_Bkg"))), 
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE),
      Normalization(normDSTot, RooAbsReal::NumEvent),
      LineColor(kGreen+3), LineStyle(1), Precision(1e-4), NumCPU(nCPU)
      );
  ws->data("dsToFit")->plotOn(myPlot2_H,Name("data_Mass"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_H,Name("Mass_Tot"),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE),
      Normalization(normDSTot, RooAbsReal::NumEvent),
      NumCPU(nCPU), LineColor(kBlack)
      );
  //ws->saveSnapshot("pdfCTAUMASS_Tot_parFit",*newpars,kTRUE);
  TH1* h = ws->data("dsToFit")->createHistogram("hist2", *ws->var("mass"), Binning(myPlot_H->GetNbinsX(), myPlot_H->GetXaxis()->GetXmin(),myPlot_H->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  // Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBiAbove(0.0)) );
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));
  double Ydown;
  double Yup;
  Ydown = YMin/(TMath::Power((YMax/YMin), (0.1/(1.0-0.1-0.4))));
  Yup = YMax*TMath::Power((YMax/YMin), (0.4/(1.0-0.1-0.4)));
  myPlot2_H->GetYaxis()->SetRangeUser(Ydown,Yup);
  myPlot2_H->GetXaxis()->SetRangeUser(2.6, 3.5);
  myPlot2_H->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-} (GeV/c^{2})}");
  myPlot2_H->SetFillStyle(4000);
  myPlot2_H->GetXaxis()->SetLabelSize(0);
  myPlot2_H->GetXaxis()->SetTitleSize(0);
  myPlot2_H->Draw();
  TLegend* leg_H = new TLegend(text_x+0.25,text_y+0.03,text_x+0.38,text_y-0.17); leg_H->SetTextSize(text_size);
  leg_H->SetTextFont(43);
  leg_H->SetBorderSize(0);
  leg_H->AddEntry(myPlot2_H->findObject("data_Mass"),"Data","pe");
  leg_H->AddEntry(myPlot2_H->findObject("Mass_Tot"),"Total fit","fl");
  leg_H->AddEntry(myPlot2_H->findObject("Mass_Bkg"),"Background","fl");
  leg_H->AddEntry(myPlot2_H->findObject("Mass_JpsiPR"),"J/#psi Prompt","l");
  leg_H->AddEntry(myPlot2_H->findObject("Mass_JpsiNoPR"),"J/#psi Non-Prompt","l");
  leg_H->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);
  drawText(Form("N_{J/#psi} = %.f #pm %.f",ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError()),text_x+0.5,text_y-y_diff,text_color,text_size);
  drawText(Form("N_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(), ws->var("N_Bkg")->getError() ),text_x+0.5,text_y-y_diff*2,text_color,text_size);
  drawText(Form("b_{J/#psi} = %.5f #pm %.5f",ws->var("b_Jpsi")->getVal(),ws->var("b_Jpsi")->getError()),text_x+0.5,text_y-y_diff*3,text_color,text_size);

  TPad *pad_H_2 = new TPad("pad_H_2", "pad_H_2", 0, 0.006, 0.98, 0.227);
  RooPlot* frameTMP_H = (RooPlot*)myPlot2_H->Clone("TMP_H");
  RooHist* hpull_H;
  pullDist(ws, pad_H_2, c_H, frameTMP_H, hpull_H, "data_Mass", "Mass_Tot", "mass", nMassBin, massLow, massHigh, "m_{#mu^{+}#mu^{-} (GeV/c^{2})}");
  pad_H_2->Update();
  printChi2(ws, pad_H_2, frameTMP_H, fitResult, "mass", "data_Mass", "Mass_Tot", nMassBin, false);

  cout << "############################################################################" << endl;
  const int NCUT = 2001;

  double IntegralCut[NCUT]={-1.00000};
  for(int i=1; i<NCUT+1; i++){
    IntegralCut[i]=IntegralCut[i-1]+0.001; //cout<<i<<endl;
  }

  const int Nctauarr = sizeof(IntegralCut)/sizeof(double);
  TH1D* outh1 = new TH1D("Fraction1","Fraction (#font[12]{l}_{J/#psi}<Cut);#font[12]{l}_{J/#psi};nonprompt J/#psi Fraction",Nctauarr-1,IntegralCut);
  TH1D* outh2 = new TH1D("Fraction2","Fraction (#font[12]{l}_{J/#psi}>Cut);#font[12]{l}_{J/#psi};nonprompt J/#psi Fraction",Nctauarr-1,IntegralCut);
  TH1D* outh3 = new TH1D("Fraction3","Fraction (Cut1<#font[12]{l}_{J/#psi}<Cut2);#font[12]{l}_{J/#psi};nonprompt J/#psi Fraction",Nctauarr-1,IntegralCut);
  TH1D* outh4 = new TH1D("Fraction4","Fraction Inclusive", 1,0,1);
  double temp1;
  double temperr1;
  double fraction=ws->var("b_Jpsi")->getVal();
  outh4->SetBinContent(0,fraction);
  outh4->SetBinError(0,ws->var("b_Jpsi")->getError());

  ofstream of1(Form("figs/2DFit_%s/Final/2DFit_Ctau_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_Left.txt", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  ofstream of2(Form("figs/2DFit_%s/Final/2DFit_Ctau_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_Righ.txt", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  ofstream of3(Form("figs/2DFit_%s/Final/2DFit_Ctau_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_Cent.txt", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));

  cout<<"test"<<endl;
  cout<<"total: "<<ws->var("N_Jpsi")->getVal()<<", 1/3: "<<ws->var("N_Jpsi")->getVal()/3<<endl;
  RooAbsPdf* test11 = (RooAbsPdf*) ws->pdf("pdfCTAUMASS_JpsiPR");
  RooAbsPdf* test22 = (RooAbsPdf*) ws->pdf("pdfCTAUMASS_JpsiNoPR");
  RooAbsReal* PR_tot = test11->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range("ctauRange"));
  RooAbsReal* NP_tot = test22->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range("ctauRange"));
  double CutBin=(PR_tot->getVal()+NP_tot->getVal())*0.35;
  double CutBin2=(PR_tot->getVal()+NP_tot->getVal())*0.35;
  cout<<"total: "<<PR_tot->getVal()+NP_tot->getVal()<<", 1/3: "<<(PR_tot->getVal()+NP_tot->getVal())/3<<endl;

  double CutLeft;
  double residualPR;
  double residualNP;
  if(ptLow==12&&ptHigh==50){residualPR=0.03;residualNP=0.05;}
  else residualPR=0.05;residualNP=0.045;
  for(int i = 0; i<NCUT; i++){
    ws->var("ctau3D")->setRange("ctauLow", ctauLow, IntegralCut[i]);
    RooAbsPdf* Low1 = (RooAbsPdf*) ws->pdf("pdfCTAUMASS_JpsiPR");
    RooAbsPdf* Low2 = (RooAbsPdf*) ws->pdf("pdfCTAUMASS_JpsiNoPR");
    RooAbsReal* PR_Integral = Low1->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range("ctauLow"));
    RooAbsReal* NP_Integral = Low2->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range("ctauLow"));
    //of1 << "[INFO] ctau Cut: "<<IntegralCut[i]<<endl;
    //of1 << "[INFO] PR Integral: "<< PR_Integral->getVal() <<endl;
    //of1 << "[INFO] NP Integral: "<< NP_Integral->getVal() <<endl;
    //of1 << "[INFO] BJpsi Fraction: "<<NP_Integral->getVal()*ws->var("b_Jpsi")->getVal()/(ws->var("b_Jpsi")->getVal()*NP_Integral->getVal()+(1-ws->var("b_Jpsi")->getVal())*PR_Integral->getVal())<<endl;
    //of1 << "[INFO] Integral: "<<NP_Integral->getVal()+PR_Integral->getVal()<<"\n"<<endl;
    //if(PR_Integral->getVal()+NP_Integral->getVal()>=CutBin) {
    //if(PR_Integral->getVal()>=0.40) {
    if(NP_Integral->getVal()>=residualNP) {
    //if(NP_Integral->getVal()*ws->var("b_Jpsi")->getVal()/(ws->var("b_Jpsi")->getVal()*NP_Integral->getVal()+(1-ws->var("b_Jpsi")->getVal())*PR_Integral->getVal())<=0.1) {
      temp1=NP_Integral->getVal()*fraction/(fraction*NP_Integral->getVal()+(1-fraction)*PR_Integral->getVal());
      outh1->SetBinContent(i,temp1);
      cout<<"Amount: "<<CutBin<<"="<<PR_Integral->getVal()+NP_Integral->getVal()<<endl;
      CutLeft=IntegralCut[i]; cout<<"CUT LOW: "<<CutLeft<<endl;
      cout<<"BJpsi Fraction: "<<NP_Integral->getVal()*ws->var("b_Jpsi")->getVal()/(ws->var("b_Jpsi")->getVal()*NP_Integral->getVal()+(1-ws->var("b_Jpsi")->getVal())*PR_Integral->getVal())<<"\n"<<endl;

      of1<<"Amount: "<<CutBin<<"="<<PR_Integral->getVal()+NP_Integral->getVal()<<endl;
      CutLeft=IntegralCut[i]; of1<<"CUT LOW: "<<CutLeft<<endl;
      of1<<"BJpsi Fraction: "<<NP_Integral->getVal()*ws->var("b_Jpsi")->getVal()/(ws->var("b_Jpsi")->getVal()*NP_Integral->getVal()+(1-ws->var("b_Jpsi")->getVal())*PR_Integral->getVal())<<"\n"<<endl;
      break;}
  }
  
  double CutRight;

  for(int i = 0; i<NCUT; i++){
    ws->var("ctau3D")->setRange("ctauHigh", IntegralCut[i], ctauHigh);
    RooAbsPdf* High1 = (RooAbsPdf*) ws->pdf("pdfCTAUMASS_JpsiPR");
    RooAbsPdf* High2 = (RooAbsPdf*) ws->pdf("pdfCTAUMASS_JpsiNoPR");
    RooAbsReal* PR_Integral = High1->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range("ctauHigh"));
    RooAbsReal* NP_Integral = High2->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range("ctauHigh"));
    //of2 << "[INFO] ctau Cut: "<<IntegralCut[i]<<endl;
    //of2 << "[INFO] PR Integral: "<< PR_Integral->getVal() <<endl;
    //of2 << "[INFO] NP Integral: "<< NP_Integral->getVal() <<endl;
    //of2 << "[INFO] BJpsi Fraction: "<<NP_Integral->getVal()*ws->var("b_Jpsi")->getVal()/(ws->var("b_Jpsi")->getVal()*NP_Integral->getVal()+(1-ws->var("b_Jpsi")->getVal())*PR_Integral->getVal())<<endl;
    //of2 << "[INFO] Integral: "<<NP_Integral->getVal()+PR_Integral->getVal()<<"\n"<<endl;
    //if(PR_Integral->getVal()+NP_Integral->getVal()<=CutBin2) {
    //if(PR_Integral->getVal()<=0.07) {
    if(PR_Integral->getVal()<=residualPR) {
    //if(NP_Integral->getVal()*ws->var("b_Jpsi")->getVal()/(ws->var("b_Jpsi")->getVal()*NP_Integral->getVal()+(1-ws->var("b_Jpsi")->getVal())*PR_Integral->getVal())>=0.8) {
      temp1=NP_Integral->getVal()*fraction/(fraction*NP_Integral->getVal()+(1-fraction)*PR_Integral->getVal());
      outh2->SetBinContent(i,temp1);
      cout<<"Amount: "<<CutBin<<"="<<PR_Integral->getVal()+NP_Integral->getVal()<<endl;
      CutRight=IntegralCut[i]; cout<<"CUT HIGH: "<<CutRight<<endl; 
      cout<<"BJpsi Fraction: "<<NP_Integral->getVal()*ws->var("b_Jpsi")->getVal()/(ws->var("b_Jpsi")->getVal()*NP_Integral->getVal()+(1-ws->var("b_Jpsi")->getVal())*PR_Integral->getVal())<<"\n"<<endl;

      of2<<"Amount: "<<CutBin<<"="<<PR_Integral->getVal()+NP_Integral->getVal()<<endl;
      CutRight=IntegralCut[i]; of2<<"CUT HIGH: "<<CutRight<<endl; 
      of2<<"BJpsi Fraction: "<<NP_Integral->getVal()*ws->var("b_Jpsi")->getVal()/(ws->var("b_Jpsi")->getVal()*NP_Integral->getVal()+(1-ws->var("b_Jpsi")->getVal())*PR_Integral->getVal())<<"\n"<<endl;
      break;}
  }

  cout << "Find cut 3\n" << endl;
  ws->var("ctau3D")->setRange("ctauFix", CutLeft, CutRight); // Input Left, Right l_jpsi
  //ws->var("ctau3D")->setRange("ctauFix", CutLeft, 0.104); // Input Left, Right l_jpsi
  RooAbsPdf* test1 = (RooAbsPdf*) ws->pdf("pdfCTAUMASS_JpsiPR");
  RooAbsPdf* test2 = (RooAbsPdf*) ws->pdf("pdfCTAUMASS_JpsiNoPR");
  RooAbsReal* PR_Integral = test1->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range("ctauFix"));
  RooAbsReal* NP_Integral = test2->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range("ctauFix"));
  cout<<"Amount: "<<CutBin<<"="<<PR_Integral->getVal()+NP_Integral->getVal()<<endl;
  cout<<"CUT LOW: "<<CutLeft<<", CUT HIGH:"<<CutRight<<endl;
  cout << "[INFO] Center BJpsi Fraction: "<<NP_Integral->getVal()*ws->var("b_Jpsi")->getVal()/(ws->var("b_Jpsi")->getVal()*NP_Integral->getVal()+(1-ws->var("b_Jpsi")->getVal())*PR_Integral->getVal())<<"\n"<<endl;

  temp1=NP_Integral->getVal()*fraction/(fraction*NP_Integral->getVal()+(1-fraction)*PR_Integral->getVal());
  outh3->SetBinContent(0,temp1);
  of3<<"Amount: "<<CutBin<<"="<<PR_Integral->getVal()+NP_Integral->getVal()<<endl;
  of3<<"CUT LOW: "<<CutLeft<<", CUT HIGH: "<<CutRight<<endl;
  of3<<"BJpsi Fraction: "<<NP_Integral->getVal()*ws->var("b_Jpsi")->getVal()/(ws->var("b_Jpsi")->getVal()*NP_Integral->getVal()+(1-ws->var("b_Jpsi")->getVal())*PR_Integral->getVal())<<"\n"<<endl;

  c_G->Update();
  //c_G->SaveAs(Form("figs/2DFit_%s/Final/2DFit_Ctau_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_dh_default.png", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  c_G->SaveAs(Form("figs/2DFit_%s/Final/2DFit_Ctau_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  //c_G->SaveAs(Form("figs/2DFit_%s/Final/2DFit_Ctau_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_dh_lamda_dss_fix_bfrc0.30_bin2.png", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  //c_G->SaveAs(Form("figs/2DFit_%s/Final/2DFit_Ctau_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_dh_lamda_dss_fix_bfrc0.30_bin2.png", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  c_H->Update();
  //c_H->SaveAs(Form("figs/2DFit_%s/Final/2DFit_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_dh_default.png", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  c_H->SaveAs(Form("figs/2DFit_%s/Final/2DFit_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  //c_H->SaveAs(Form("figs/2DFit_%s/Final/2DFit_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_dh_lamda_dss_fix_bfrc0.30.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  //c_H->SaveAs(Form("figs/2DFit_%s/Final/2DFit_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_dh_lamda_dss_fix_bfrc0.30.png", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));

  TH1D* outh = new TH1D("2DfitResults","fit result",20,0,20);

  float temp = ws->var("b_Jpsi")->getVal();
  float temperr = ws->var("b_Jpsi")->getError();

  outh->SetBinContent(1,temp);
  outh->SetBinError(1,temperr);

  TFile *outFile = new TFile(Form("roots/2DFit_%s/Final/2DFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP),"recreate");
  //ws->Write();
  outFile->cd();
  outh->Write();
  outh1->Write();
  outh2->Write();
  outh3->Write();
  outh4->Write();
  outFile->Close();

  //fitResult->Print("v");
}
