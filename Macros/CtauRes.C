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

void CtauRes_up(
    float ptLow=3, float ptHigh=4.5,
    float yLow=1.6, float yHigh=2.4,
    int cLow=20, int cHigh=120,
    float muPtCut=0.0,
    int PR=2
    )
{
  gStyle->SetEndErrorSize(0);
  TFile* fMass; TFile* fCErr;
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

  TString bCont;
  if(PR==0) bCont="Prompt";
  else if(PR==1) bCont="NonPrompt";
  else if(PR==2) bCont="Inclusive";

  fMass = new TFile(Form("../roots/2DFit_1127/MassFitResult_%s_%s.root",bCont.Data(),kineLabel.Data()));
  fCErr = new TFile(Form("../roots/2DFit_1127/CtauErrResult_%s_%s.root",bCont.Data(),kineLabel.Data()));
  
  RooDataSet *dataw_Sig = (RooDataSet*)fCErr->Get("dataw_Sig");
  RooDataSet *datasetMass = (RooDataSet*)fMass->Get("datasetMass");
  RooHistPdf *pdfCTAUERR_Jpsi = (RooHistPdf*)fCErr->Get("pdfCTAUERR_Jpsi");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*datasetMass);
  ws->import(*dataw_Sig);
  ws->import(*pdfCTAUERR_Jpsi);
  RooDataSet *ctauResCutDS =(RooDataSet*)dataw_Sig->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")),*(ws->var("ctau3DErr"))),"ctau3DRes<0");
  ctauResCutDS->SetName("ctauResCutDS");
  ws->import(*ctauResCutDS);
  //RooArgSet *newpars = ws->pdf("pdfCTAUERR_Jpsi")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr"), *ws->var("ctau3DRes")));
  RooArgSet* params = (RooArgSet*) ws->pdf("pdfCTAUERR_Jpsi")->getParameters(RooArgSet(*ws->var("ctau3DRes"), *ws->var("ctau3DErr"), *ws->var("ctau3D")));
  cout << "####################################" << endl;
  //ws->var("mass")->setRange(massLow, massHigh);
  //ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  //ws->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);
  //ws->var("ctau3DErr")->setRange(ctauErrLow, ctauErrHigh);
  //ws->var("ctau3DErr")->setRange("ctauErrRange",ctauErrLow, ctauErrHigh);
  //ws->var("ctau3DRes")->setRange(ctauResLow, ctauResHigh);
  //ws->var("ctau3DRes")->setRange("ctauResRange", ctauResLow, ctauResHigh);
  //ws->var("mass")->Print();
  //ws->var("ctau3D")->Print();
  //ws->var("ctau3DErr")->Print();
  //***********************************************************************
  //**************************** CTAU RES FIT *****************************
  //***********************************************************************
  cout << endl << "*************** Start Res Fit *****************" << endl << endl;
  Double_t outTot = ws->data("dataw_Sig")->numEntries();
  Double_t outErr = outTot - (*ws->data("ctauResCutDS")->numEntries());
  //double entries = ws->data("ctauResCutDS")->sumEntries();
  //double entries = ws->var("N_Jpsi")->getVal();
  //cout<<"[Info] #J/psi: "<<entries<<endl;
  //ws->factory(Form("N_Jpsi[%.12f,%.12f,%.12f]", entries, entries, entries*2.0));
  cout << "HERE" << endl;
  // create the variables for this model
  //ws->var("ctau3DRes")->setRange("ctauResRange", ctauResLow, ctauResHigh);
  int nGauss = 3;
  ws->factory("ctauRes_mean[0.0]");
  ws->factory("ctau1_CtauRes[0., -0.1, 0.1]");  ws->factory("s1_CtauRes[1., 0., 10.]");
  ws->factory("ctau2_CtauRes[0., -0.1, 0.1]");  ws->factory("s2_CtauRes[2., 0., 10.]");
  ws->factory("ctau3_CtauRes[0., -0.1, 0.1]");  ws->factory("s3_CtauRes[3,  0., 10.]");
  ws->factory("ctau4_CtauRes[0., -0.1, 0.1]");  ws->factory("s4_CtauRes[5.37, 0., 10.]");
  ws->factory("f2_CtauRes[0.5, 0., 1.]"); ws->factory("f_CtauRes[0.5, 0., 1.]"); ws->factory("f3_CtauRes[0.5, 0., 1.]");
  // create the three PDFs
  TString varName="ctau3DRes";
  ws->factory(Form("GaussModel::%s(%s, %s, %s)", "GaussModel1_ctauRes", varName.Data(),
        "ctau1_CtauRes", //"ctau1_CtauRes",
        "s1_CtauRes"
        //"One", //meanSF
        //"SF_sigma"//sigmaSF (usePerEventError?"ctauErr":"One")
        ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s)", "GaussModel2_ctauRes", varName.Data(),
        "ctau2_CtauRes", //"ctau2_CtauRes",
        "s2_CtauRes"
        //"One", //meanSF
        //"SF_sigma"//sigmaSF (usePerEventError?"ctauErr":"One")
        ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s)", "GaussModel3_ctauRes", varName.Data(),
        "ctau3_CtauRes", //"ctau3_CtauRes",
        "s3_CtauRes"
        //"One", //meanSF
        //"SF_sigma"//sigmaSF (usePerEventError?"ctauErr":"One")
        ));
  ws->factory(Form("GaussModel::%s(%s, %s, %s)", "GaussModel4_ctauRes", varName.Data(),
        "ctau4_CtauRes", //"ctau3_CtauRes",
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

  RooAddPdf* GaussModel_Tot = new RooAddPdf("GaussModel_Tot");
  RooAbsPdf *ctauResModel = ctauResModel = new RooAddPdf("GaussModel_Tot", "GaussModel_Tot", *ws->pdf("GaussModelCOND_ctauRes"), *ws->var("N_Jpsi"));
  ws->import(*ctauResModel);

  TCanvas* c_C =  new TCanvas("canvas_C","My plots",1108,4,550,520);
  c_C->cd();
  TPad *pad_C_1 = new TPad("pad_C_1", "pad_C_1", 0, 0.16, 0.98, 1.0);
  pad_C_1->SetTicks(1,1);
  pad_C_1->Draw(); pad_C_1->cd();
  RooPlot* myPlot_C = ws->var("ctau3DRes")->frame(nCtauResBins); // bins
  myPlot_C->SetTitle("");

  pad_C_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot2_C = (RooPlot*)myPlot_C->Clone();
  bool isWeighted = ws->data("ctauResCutDS")->isWeighted();
  RooFitResult* fitCtauRes = ws->pdf("GaussModel_Tot")->fitTo(*ctauResCutDS, Save(), Range(ctauResLow,0), SumW2Error(isWeighted), Extended(kTRUE), NumCPU(10), PrintLevel(-1));
  //RooFitResult* fitCtauRes = ws->pdf("GaussModel_Tot")->fitTo(*ctauResCutDS, Save(), SumW2Error(isWeighted), Extended(kTRUE), NumCPU(10), PrintLevel(-1));
  ws->import(*fitCtauRes);
  //setFixedVarsToContantVars(ws);
  ws->data("ctauResCutDS")->plotOn(myPlot2_C, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0),
      MarkerSize(.7), LineColor(kBlack), MarkerColor(kBlack));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_ctauRes"), Precision(1e-4), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), LineColor(kBlack));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm1"), Precision(1e-4), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel1_ctauRes")), LineColor(kGreen+2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm2"), Precision(1e-4), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel2_ctauRes")), LineColor(kRed+2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm3"), Precision(1e-4), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel3_ctauRes")), LineColor(kBlue+2));
  if(nGauss==4){ws->pdf("GaussModel_Tot")->plotOn(myPlot2_C,Name("modelHist_gm4"), Precision(1e-4), Range("ctauResRange"),
      Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel4_ctauRes")), LineColor(kMagenta+2));}
  ws->data("ctauResCutDS")->plotOn(myPlot2_C, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0),
      MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2), Precision(1e-4));
  
  TH1* h = ws->data("ctauResCutDS")->createHistogram("hist", *ws->var("ctau3DRes"), Binning(myPlot_C->GetNbinsX(),myPlot_C->GetXaxis()->GetXmin(),myPlot_C->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));
  Double_t Yup(0.),Ydown(0.);
  Yup = YMax*TMath::Power((YMax/YMin), (0.5/(1.0-0.5-0.2)));
  Ydown = YMin/(TMath::Power((YMax/YMin), (0.2/(1.0-0.5-0.2))));
  myPlot2_C->GetYaxis()->SetRangeUser(Ydown,Yup);
  cout<<Ydown<<" to "<<Yup<<endl;
  cout<<"###"<<endl;
  cout<<ws->var("ctau3DErr")->getMin()<<", "<<ws->var("ctau3DErr")->getMax()<<endl;
  vector<double> rangeErr; rangeErr.push_back(ws->var("ctau3DErr")->getMin()); rangeErr.push_back(ws->var("ctau3DErr")->getMax());
  if (outErr>0.0) {
    TLine   *minline = new TLine(rangeErr[0], 0.0, rangeErr[0], Ydown*TMath::Power((Yup/Ydown),0.4));
    minline->SetLineStyle(2);
    minline->SetLineColor(2);
    minline->SetLineWidth(3);
    myPlot2_C->addObject(minline);
    TLine   *maxline = new TLine(rangeErr[1], 0.0, rangeErr[1], Ydown*TMath::Power((Yup/Ydown),0.4));
    maxline->SetLineStyle(2);
    maxline->SetLineColor(1);
    maxline->SetLineWidth(3);
    myPlot2_C->addObject(maxline);}
  myPlot2_C->GetXaxis()->CenterTitle();
  myPlot2_C->GetXaxis()->SetTitle("#frac{l_{J/#psi}}{#sigma_{l_{J/#psi}}}");
  myPlot2_C->SetFillStyle(4000);
  myPlot2_C->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_C->GetXaxis()->SetLabelSize(0);
  myPlot2_C->GetXaxis()->SetTitleSize(0);
  myPlot2_C->Draw();

  TLegend* leg_C = new TLegend(text_x+0.29,text_y+0.03,text_x+0.39,text_y-0.17); leg_C->SetTextSize(text_size);
  leg_C->SetTextFont(43);
  leg_C->SetBorderSize(0);
  leg_C->AddEntry(myPlot2_C->findObject("dataHist_ctauRes"),"Data","pe");
  leg_C->AddEntry(myPlot2_C->findObject("modelHist_ctauRes"),"Total PDF","l");
  leg_C->AddEntry(myPlot2_C->findObject("modelHist_gm1"),"Gauss 1","l");
  leg_C->AddEntry(myPlot2_C->findObject("modelHist_gm2"),"Gauss 2","l");
  leg_C->AddEntry(myPlot2_C->findObject("modelHist_gm3"),"Gauss 3","l");
  if(nGauss==4)leg_C->AddEntry(myPlot2_C->findObject("modelHist_gm4"),"Gauss 4","l");
  leg_C->Draw("same");
  cout<<"s2/s1: "<<ws->var("s2_CtauRes")->getVal()/ws->var("s1_CtauRes")->getVal()<<endl;
  cout<<"s3/s2: "<<ws->var("s3_CtauRes")->getVal()/ws->var("s2_CtauRes")->getVal()<<endl;
  cout<<"s1: "<<ws->var("s1_CtauRes")->getVal()<<endl;
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
  if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
  else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
  drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);
  float dy = 0;
  drawText(Form("Excl: (%.4f%%) %.2f evts", (outErr*100.0/outTot), outErr),text_x,text_y-y_diff*3,text_color,text_size); dy+=1.5*0.045;

  drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError()),text_x+0.5,text_y,text_color,text_size);
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

  pullFrame_C->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#font[12]{#sigma}_{J/#psi}}");
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
  RooArgSet* fitargs = new RooArgSet();
  fitargs->add(fitCtauRes->floatParsFinal());
  RooDataSet *datasetRes = new RooDataSet("datasetRes","dataset with Resolution Fit result", *fitargs);

  c_C->Update();
  c_C->SaveAs(Form("../figs/2DFit_1127/CtauRes_%s_%s.pdf",bCont.Data(),kineLabel.Data()));

  TFile *outFile = new TFile(Form("../roots/2DFit_1127/CtauResResult_%s_%s.root",bCont.Data(),kineLabel.Data()),"recreate");
  ctauResModel->Write();
  //GaussModel_Tot->Write();
  //ctauResCutDS->Write();
  //datasetRes->Write();
  fitCtauRes->Write();
  outFile->Close();
}
