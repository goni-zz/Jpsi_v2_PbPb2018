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

void CtauBkg_up(
    float ptLow=12, float ptHigh=50,
    float yLow=1.6, float yHigh=2.4,
    int cLow=20, int cHigh=120,
    float muPtCut=0.0,
    int PR=2, //0=PR, 1=NP, 2=Inc.
    float ctau3DErrCut=0.14
    )
{

  gStyle->SetEndErrorSize(0);
  gSystem->mkdir("../roots/2DFit_1127/");
  gSystem->mkdir("../figs/2DFit_1127/");

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

  TFile* f1; TFile* fMass; TFile* fCErr; TFile* fCRes;
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

  f1 = new TFile("../skimmedFiles/OniaRooDataSet_isMC0_JPsi_w_Effw0_Accw0_PtW1_TnP1_20201127.root");
  fMass = new TFile(Form("../roots/2DFit_1127/MassFitResult_%s_%s.root",bCont.Data(),kineLabel.Data()));
  fCErr = new TFile(Form("../roots/2DFit_1127/CtauErrResult_%s_%s.root",bCont.Data(),kineLabel.Data()));
  fCRes = new TFile(Form("../roots/2DFit_1127/CtauResResult_%s_%s.root",bCont.Data(),kineLabel.Data()));

  RooDataSet *datasetMass = (RooDataSet*)fMass->Get("datasetMass");
  RooAddPdf* pdfMASS_Tot = (RooAddPdf*)fMass->Get("pdfMASS_Tot");
  RooDataSet *dataw_Bkg = (RooDataSet*)fCErr->Get("dataw_Bkg");
  RooAddPdf* GaussModel_Tot = (RooAddPdf*)fCRes->Get("GaussModel_Tot");

  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*datasetMass);
  ws->import(*pdfMASS_Tot);
  ws->import(*dataw_Bkg);
  ws->import(*GaussModel_Tot);
  cout<<"####################################" <<endl;
  cout<<"pt: "<<ptLow<<"-"<<ptHigh<<", y: "<<yLow<<"-"<<yHigh<<", Cent: "<<cLow/2<<"-"<<cHigh/2<<"%"<<endl;

  cout <<"******** New Combined Dataset ***********" <<endl;
  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  ws->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);
  ws->var("ctau3DErr")->setRange(ctauErrLow, ctauErrHigh);
  ws->var("ctau3DErr")->setRange("ctauErrRange",ctauErrLow, ctauErrHigh);
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
  if(ptLow==3&&ptHigh==4.5){
  ws->factory("b_Bkg[0.1, 1e-6, 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.5, 1e-6, 1.]");
  ws->factory("fDLIV[0.5, 1e-6, 1.]");
  ws->factory("lambdaDDS_Bkg[0.5, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg[ 0., 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg[0.5, 1e-6, 1.]");}
  if(ptLow==4.5&&ptHigh==6.5){
  ws->factory("b_Bkg[1e-6, 1e-6, 1.]");//NP fraction for bkg
  ws->factory("fDFSS[1e-6, 1e-6, 1.]");
  ws->factory("fDLIV[1e-6, 1e-6, 1.]");
  ws->factory("lambdaDDS_Bkg[1e-6, 0.5, 1.]");
  ws->factory("lambdaDF_Bkg[ 1e-6.,0.03, 1.]");
  ws->factory("lambdaDSS_Bkg[1e-6, 0.06, 1.]");}
  if(ptLow==6.5&&ptHigh==9){
  ws->factory("b_Bkg[0.6, 1e-6, 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.5, 1e-6, 1.]");
  ws->factory("fDLIV[0.5, 1e-6, 1.]");
  ws->factory("lambdaDDS_Bkg[0.2, 0.0, 0.4]");
  ws->factory("lambdaDF_Bkg[0.3, 0.1, 1.]");
  ws->factory("lambdaDSS_Bkg[0.5, 1e-6, 1.]");}
  if(ptLow==9&&ptHigh==12){
  ws->factory("b_Bkg[0.5, 1e-6, 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.5, 1e-6, 1.]");
  ws->factory("fDLIV[0.5, 1e-6, 1.]");
  ws->factory("lambdaDDS_Bkg[0.2, 1e-6, 0.4]");
  ws->factory("lambdaDF_Bkg[0.3, 0.1, 1.]");
  ws->factory("lambdaDSS_Bkg[0.5, 1e-6, 1.]");}
  if(ptLow==12&&ptHigh==50){
  ws->factory("b_Bkg[0.1, 0., 1.]");//NP fraction for bkg
  ws->factory("fDFSS[0.5, 1e-6, 1.]");
  ws->factory("fDLIV[0.5, 1e-6, 1.]");
  ws->factory("lambdaDDS_Bkg[0.2, 1e-6, 1.]");
  ws->factory("lambdaDF_Bkg[0.3, 1e-6, 1.]");
  ws->factory("lambdaDSS_Bkg[0.5, 1e-6, 1.]");}
  //parameters fixed by Resolution model
  ws->var("ctau1_CtauRes")->setConstant(kTRUE); ws->var("s1_CtauRes")->setConstant(kTRUE);
  ws->var("ctau2_CtauRes")->setConstant(kTRUE);	ws->var("s2_CtauRes")->setConstant(kTRUE);
  ws->var("ctau3_CtauRes")->setConstant(kTRUE);	ws->var("s3_CtauRes")->setConstant(kTRUE);
  ws->var("f2_CtauRes")->setConstant(kTRUE);	ws->var("f_CtauRes")->setConstant(kTRUE);
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
  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes3", "ctau3D",
        "ctau3_CtauRes", //"ctau3_CtauRes",
        "s3_CtauRes",
        "zeroMean",
        "ctau3DErr"
        ));
  ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "ctauRes32", "ctauRes3", "ctauRes2", "f2_CtauRes"));
  ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "pdfCTAURES", "ctauRes32", "ctauRes1", "f_CtauRes"));
  //make 3 exp
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS", "ctau3D", "lambdaDSS_Bkg", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF", "ctau3D", "lambdaDF_Bkg", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", "pdfCTAUDDS", "ctau3D", "lambdaDDS_Bkg", "pdfCTAURES"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU1", "fDFSS", "pdfCTAUDSS", "pdfCTAUDF"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_BkgNoPR", "fDLIV", "pdfCTAU1", "pdfCTAUDDS"));//NP
  ws->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_BkgPR","pdfCTAURES"));//PR
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_Bkg", "b_Bkg", "pdfCTAUCOND_BkgNoPR", "pdfCTAUCOND_BkgPR"));
  ws->factory(Form("RooExtendPdf::%s(%s,%s)", "pdfTot_Bkg", "pdfCTAUCOND_Bkg", "N_Bkg"));//N_Bkg is number of bkg from dataw_Bkg
  RooAddPdf* pdfTot_Bkg = new RooAddPdf("pdfTot_Bkg");
  //RooAbsPdf* ctauBkgModel = ctauBkgModel = new RooAddPdf("pdfTot_Bkg","pdfTot_Bkg",*ws->pdf("pdfCTAUCOND_Bkg"), *ws->var("N_Bkg"));
  pdfTot_Bkg = new RooAddPdf("pdfTot_Bkg","",RooArgList(*ws->pdf("pdfCTAUCOND_Bkg")),RooArgList(*ws->var("N_Bkg")));
  //ws->import(*pdfCTAUCOND_Bkg);

  TCanvas* c_E =  new TCanvas("canvas_E","My plots",554,565,550,520);
  c_E->cd();
  TPad *pad_E_1 = new TPad("pad_E_1", "pad_E_1", 0, 0.16, 0.98, 1.0);
  pad_E_1->SetTicks(1,1);
  pad_E_1->Draw(); pad_E_1->cd();
  RooPlot* myPlot_E = ws->var("ctau3D")->frame(Bins(nCtauBins), Range(ctauLow, ctauHigh)); // bins
  myPlot_E->SetTitle("");

  pad_E_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot2_E = (RooPlot*)myPlot_E->Clone();
  //RooFitResult* fitCtauBkg = ws->pdf("pdfTot_Bkg")->fitTo(*dataw_Bkg, Save(), Range("ctauRange"), Extended(kTRUE), NumCPU(nCPU), PrintLevel(-1));
  RooFitResult* fitCtauBkg = ws->pdf("pdfTot_Bkg")->fitTo(*dataw_Bkg, Save(), Range("ctauRange"), Extended(kTRUE), NumCPU(nCPU), PrintLevel(-1));
  ws->data("dataw_Bkg")->plotOn(myPlot2_E,Name("data_ctauBkg"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kBlue+2), MarkerColor(kBlue+2));
  ws->pdf("pdfTot_Bkg")->plotOn(myPlot2_E,Name("ctauBkg_Tot"),
      Normalization(ws->data("dataw_Bkg")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE),
      FillStyle(1001), FillColor(kAzure-9), LineColor(kBlack), LineWidth(3), DrawOption("LCF"), Precision(1e-4));
  ws->pdf("pdfTot_Bkg")->plotOn(myPlot2_E,Name("pdfCTAUCOND_BkgPR"),
      Normalization(ws->data("dataw_Bkg")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
      ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE), Components("pdfCTAUCOND_BkgPR"),
      LineWidth(2), LineColor(kRed+2));
  ws->data("dataw_Bkg")->plotOn(myPlot2_E,Name("data_ctauBkg"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kBlue+2), MarkerColor(kBlue+2));
  myPlot2_E->GetYaxis()->SetRangeUser(10e-2, 10e7);
  TH1* h = ws->data("dataw_Bkg")->createHistogram("hist", *ws->var("ctau3D"), Binning(myPlot_E->GetNbinsX(),myPlot_E->GetXaxis()->GetXmin(),myPlot_E->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));
  Double_t Yup(0.),Ydown(0.);
  Yup = YMax*TMath::Power((YMax/0.1), 0.5);
  Ydown = 0.1;
  myPlot2_E->GetYaxis()->SetRangeUser(Ydown,Yup);
  myPlot2_E->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  myPlot2_E->SetFillStyle(4000);
  myPlot2_E->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_E->GetXaxis()->SetLabelSize(0);
  myPlot2_E->GetXaxis()->SetTitleSize(0);
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
  drawText(Form("#lambdaDF_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDF_Bkg")->getVal(), ws->var("lambdaDF_Bkg")->getError() ),text_x+0.5,text_y-y_diff*5,text_color,text_size);
  drawText(Form("#lambdaDSS_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDSS_Bkg")->getVal(), ws->var("lambdaDSS_Bkg")->getError() ),text_x+0.5,text_y-y_diff*6,text_color,text_size);
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
  c_E->SaveAs(Form("../figs/2DFit_1127/Bkg_%s_%s.pdf",bCont.Data(),kineLabel.Data()));
  RooArgSet* fitargs = new RooArgSet();
  fitargs->add(fitCtauBkg->floatParsFinal());
  RooDataSet *datasetCBkg = new RooDataSet("datasetCBkg","dataset with Ctau Background Fit result", *fitargs);
  RooWorkspace *wscbkg = new RooWorkspace("workspaceCBkg");
  wscbkg->import(*fitCtauBkg);
  //wscbkg->import(*pdfCTAUCOND_Bkg);

  //	ws->Print();

  TFile *outFile = new TFile(Form("../roots/2DFit_1127/CtauBkgResult_%s_%s.root",bCont.Data(),kineLabel.Data()),"recreate");
  fitCtauBkg->Write();
  pdfTot_Bkg->Write();
  //pdfCTAUCOND_Bkg->Write();
  //pdfTot_Bkg->Write();
  datasetCBkg->Write();
  wscbkg->Write();
  outFile->Close();
}
