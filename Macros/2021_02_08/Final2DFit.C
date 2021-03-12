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

void Final2DFit(
    float ptLow=3, float ptHigh=4.5,
    float yLow=1.6, float yHigh=2.4,
    int cLow=20, int cHigh=120,
    int PR=2 //0=PR, 1=NP, 2=Inc.
    )
{
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir("../../roots/2DFit_20210310/", kTRUE);
  gSystem->mkdir("../../figs/2DFit_20210310/", kTRUE);

  TString bCont;
  if(PR==0) {bCont="Prompt";
    cout<<"#### Prompt ENHANCED"<<endl;}
  else if(PR==1) bCont="NonPrompt";
  else if(PR==2) {bCont="Inclusive";
    cout<<"#### Inclusive Jpsi"<<endl;}
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

  f1 = new TFile("../../skimmedFiles/OniaRooDataSet_isMC0_JPsi_w_Effw0_Accw0_PtW1_TnP1_20210111.root");
  //f1 = new TFile("../skimmedFiles/OniaRooDataSet_isMC0_JPsi_w_Effw0_Accw0_PtW1_TnP1_202020210310.root");
  fMass = new TFile(Form("../../roots/2DFit_20210310/Mass/MassFitResult_%s_%s.root",bCont.Data(),kineLabel.Data()));
  fCErr = new TFile(Form("../../roots/2DFit_20210310/CtauErr/CtauErrResult_%s_%s.root",bCont.Data(),kineLabel.Data()));
  fCRes = new TFile(Form("../../roots/2DFit_20210310/CtauRes/CtauResResult_%s_%s.root",bCont.Data(),kineLabel.Data()));
  //fCBkg = new TFile(Form("../../roots/2DFit_20210310/CtauBkg/CtauBkgResult_%s_%s.root",bCont.Data(),kineLabel.Data()));
  fCBkg = new TFile(Form("../roots/2DFit_20210208/CtauBkg/CtauBkgResult_%s_%s.root",bCont.Data(),kineLabel.Data()));
  fCTrue = new TFile(Form("../../roots/2DFit_20210310/CtauTrue/CtauTrueResult_%s_%s.root",bCont.Data(),kineLabel.Data()));

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
  ws->data("dataw_Bkg")->getRange(*ws->var("ctau3DErr"), ctauErrMin, ctauErrMax);
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
  RooDataSet *dsTot = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  //RooDataSet *dsToFit = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), Form("ctau3DErr>=%.6f&&ctau3DErr<=%.6f&&%s", ctauErrMin, ctauErrMax, kineCut.Data()) );
  RooDataSet* dsToFit = (RooDataSet*)dsTot->reduce(Form("ctau3DErr>=%.6f&&ctau3DErr<=%.6f",ctauErrMin, ctauErrMax))->Clone("dsTot");
  dsToFit->SetName("dsToFit");
  dsTot->SetName("dsTot");
  ws->import(*dsTot);
  ws->import(*dsToFit);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  ws->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);
  ws->var("ctau3DErr")->setRange(ctauErrMin, ctauErrMax);
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
  double fdss1 = ws->var("f")->getVal();
  //double lambda2 = ws->var("lambdaDSS2")->getVal();
  //ws->var("b_Bkg")->setConstant(kTRUE);//
  //make jpsi pdf
  ws->factory(Form("lambdaDSS_test1[%.4f, %.4f, %.4f]", lambda, 1e-8, lambda*2));
  ws->factory(Form("lambdaDSS_test2[%.4f, %.4f, %.4f]", lambda1, 1e-8, lambda1*2));
  ws->factory(Form("lambdaDSS_test3[%.4f, %.4f, %.4f]", lambda2, 1e-8, lambda2*2));
  ws->factory(Form("fDSS1[%.4f, %.4f, %.4f]", fdss, 1e-8, 1.));
  ws->factory(Form("fDSS2[%.4f, %.4f, %.4f]", fdss1, 1e-8, 1.));
  //NoPR{
  //1exp
  //ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUCOND_JpsiNoPR", "ctau3D", "lambdaDSS", "pdfCTAURES")); //NP
  //3exp
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUE_test1", "ctau3D", "lambdaDSS_test1", "pdfCTAURES")); //NP
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUE_test2", "ctau3D", "lambdaDSS_test2", "pdfCTAURES")); //NP
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUE_test3", "ctau3D", "lambdaDSS_test3", "pdfCTAURES")); //NP
  ws->factory(Form("AddModel::%s({%s , %s}, %s)", "pdfCTAUTRUE_test12", "pdfCTAUTRUE_test1", "pdfCTAUTRUE_test2", "fDSS"));
  ws->factory(Form("AddModel::%s({%s , %s}, %s)", "pdfCTAUCOND_JpsiNoPR", "pdfCTAUTRUE_test12", "pdfCTAUTRUE_test3", "fDSS1"));//}
  //PR
  ws->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_JpsiPR", "pdfCTAURES"));
  //3-4.5
  if(ptLow==3&&ptHigh==4.5){ws->factory("b_Jpsi[0.2027, 1e-6, 1.]");}//NP fraction for Sig
  else if(ptLow==6.5&&ptHigh==9){ws->factory("b_Jpsi[0.19, 0.19, .4]");}//NP fraction for Sig
  else if(ptLow==6.5&&ptHigh==8.5&&yHigh==0.6){ws->factory("b_Jpsi[0.27, 0.21, .29]");}//NP fraction for Sig
  else if(ptLow==9&&ptHigh==10){ws->factory("b_Jpsi[0.27, 0.0, 1.]");}//NP fraction for Sig
  //ws->factory("b_Jpsi[0.27, 1e-8, 1.]");//NP fraction for Sig

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

  TCanvas* c_G =  new TCanvas("canvas_G","My plots",1108,4,550,520);
  c_G->cd();
  TPad *pad_G_1 = new TPad("pad_G_1", "pad_G_1", 0, 0.16, 0.98, 1.0);
  pad_G_1->SetTicks(1,1);
  pad_G_1->Draw(); pad_G_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot_G = ws->var("ctau3D")->frame(nCtauBins); // bins
  //RooPlot* myPlot_H = ws->var("mass")->frame(Bins(nMassBin), Range(2.6,3.5)); // bins
  myPlot_G->SetTitle("");

  c_G->cd();
  c_G->SetLogy();
  pad_G_1->cd();

  double normDSTot = ws->data("dsToFit")->sumEntries()/ws->data("dsTot")->sumEntries();
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
  drawText(Form("b_{J/#psi} = %.4f #pm %.4f",ws->var("b_Jpsi")->getVal(),ws->var("b_Jpsi")->getError()),text_x+0.5,text_y+0.05-y_diff*3,text_color,text_size);

  TPad *pad_G_2 = new TPad("pad_G_2", "pad_G_2", 0, 0.006, 0.98, 0.227);
  RooPlot* frameTMP_G = (RooPlot*)myPlot2_G->Clone("TMP_G");
  RooHist* hpull_G;
  pullDist(ws, pad_G_2, c_G, frameTMP_G, hpull_G, "data_Ctau", "Ctau_Tot", "ctau3D", nCtauBins, -4, 6.0, "#font[12]{l}_{J/#psi} (mm)");
  printChi2(ws, pad_G_2, frameTMP_G, fitResult, "ctau3D", "data_Ctau", "Ctau_Tot", nCtauBins);
  pad_G_2->Update();

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
      FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), NumCPU(nCPU), LineStyle(kDashed)
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
  drawText(Form("b_{J/#psi} = %.4f #pm %.4f",ws->var("b_Jpsi")->getVal(),ws->var("b_Jpsi")->getError()),text_x+0.5,text_y-y_diff*3,text_color,text_size);

  TPad *pad_H_2 = new TPad("pad_H_2", "pad_H_2", 0, 0.006, 0.98, 0.227);
  RooPlot* frameTMP_H = (RooPlot*)myPlot2_H->Clone("TMP_H");
  RooHist* hpull_H;
  pullDist(ws, pad_H_2, c_H, frameTMP_H, hpull_H, "data_Mass", "Mass_Tot", "mass", nMassBin, massLow, massHigh, "m_{#mu^{+}#mu^{-} (GeV/c^{2})}");
  pad_H_2->Update();
  printChi2(ws, pad_H_2, frameTMP_H, fitResult, "mass", "data_Mass", "Mass_Tot", nMassBin, false);

  const int nCut = 15;
  float IntegralCut[nCut]={0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
  for(int i = 0; i<nCut; i++){
    ws->var("ctau3D")->setRange("ctauLeft", ctauLow, IntegralCut[i]);
    RooAbsPdf* test1 = (RooAbsPdf*) ws->pdf("pdfCTAUMASS_JpsiPR");
    RooAbsPdf* test2 = (RooAbsPdf*) ws->pdf("pdfCTAUMASS_JpsiNoPR");
    RooAbsReal* PR_Integral = test1->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range("ctauLeft"));
    RooAbsReal* NP_Integral = test2->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range("ctauLeft"));

    cout << "[INFO] ctau Cut: "<<IntegralCut[i]<<endl;
    cout << "[INFO] PR Integral: "<< PR_Integral->getVal() <<endl;
    cout << "[INFO] NP Integral: "<< NP_Integral->getVal() <<endl;
    cout << "[INFO] BJpsi Fraction: "<<NP_Integral->getVal()*ws->var("b_Jpsi")->getVal()/(ws->var("b_Jpsi")->getVal()*NP_Integral->getVal()+(1-ws->var("b_Jpsi")->getVal())*PR_Integral->getVal())<<"\n"<<endl;
  }

  for(int i = 0; i<nCut; i++){
    ws->var("ctau3D")->setRange("ctauRight", IntegralCut[i], ctauHigh);
    RooAbsPdf* test1 = (RooAbsPdf*) ws->pdf("pdfCTAUMASS_JpsiPR");
    RooAbsPdf* test2 = (RooAbsPdf*) ws->pdf("pdfCTAUMASS_JpsiNoPR");
    RooAbsReal* PR_Integral = test1->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range("ctauRight"));
    RooAbsReal* NP_Integral = test2->createIntegral(*ws->var("ctau3D"), NormSet(*ws->var("ctau3D")), Range("ctauRight"));

    cout << "[INFO] ctau Cut: "<<IntegralCut[i]<<endl;
    cout << "[INFO] PR Integral: "<< PR_Integral->getVal() <<endl;
    cout << "[INFO] NP Integral: "<< NP_Integral->getVal() <<endl;
    cout << "[INFO] BJpsi Fraction: "<<NP_Integral->getVal()*ws->var("b_Jpsi")->getVal()/(ws->var("b_Jpsi")->getVal()*NP_Integral->getVal()+(1-ws->var("b_Jpsi")->getVal())*PR_Integral->getVal())<<"\n"<<endl;
  }

  c_G->Update();
  c_G->SaveAs(Form("../../figs/2DFit_20210310/Ctau_%s_%s.pdf",bCont.Data(),kineLabel.Data()));
  c_H->Update();
  c_H->SaveAs(Form("../../figs/2DFit_20210310/Mass_%s_%s.pdf",bCont.Data(),kineLabel.Data()));
  TFile *outFile = new TFile(Form("../../roots/2DFit_20210310/2DFitResult_%s_%s.root",bCont.Data(),kineLabel.Data()),"recreate");
  ws->Write();
  outFile->Close();
}
