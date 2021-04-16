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
void CtauErr(
    float ptLow=6.5, float ptHigh=7.5,
    float yLow=0, float yHigh=2.4,
    int cLow=20, int cHigh=120,
    int PR=2, // 0=PR, 1=NP, 2=Inclusive
    int ordC=4
    )
{
  TString DATE="20210208";
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir(Form("../../roots/2DFit_%s/CtauErr",DATE.Data()), kTRUE);
  gSystem->mkdir(Form("../../figs/2DFit_%s/CtauErr",DATE.Data()), kTRUE);

  TString bCont;
  if(PR==0) bCont="Prompt";
  else if(PR==1) bCont="NonPrompt";
  else if(PR==2) bCont="Inclusive";

  TString fOrd;
  if (ordC==1) fOrd="1stCheby";
  else if (ordC==2) fOrd="2ndCheby";
  else if (ordC==3) fOrd="3rdCheby";
  else if (ordC==4) fOrd="4thCheby";

  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  TFile* f1; TFile* f2; TFile* f3; TFile* fMass;
  TString kineCut;
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);

  //f1 = new TFile("../skimmedFiles/OniaRooDataSet_isMC0_JPsi1SW_20201127.root");
  f1 = new TFile("../../skimmedFiles/OniaRooDataSet_isMC0_JPsi_w_Effw0_Accw0_PtW1_TnP1_20210111.root");
  //f1 = new TFile("../skimmedFiles/OniaRooDataSet_isMC0_JPsi_w_Effw0_Accw0_PtW1_TnP1_20201127.root");
  //fMass = new TFile(Form("../roots/2DFit_20210208/MassFitResult_%s_%s_%s.root",bCont.Data(),kineLabel.Data(),fOrd.Data()));
  fMass = new TFile(Form("../../roots/2DFit_20210208/Mass/MassFitResult_%s_%s.root",bCont.Data(),kineLabel.Data()));
  kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.5 && cBin>%d && cBin<%d",ptLow, ptHigh, yLow, yHigh, cLow, cHigh);

  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

  TString OS="recoQQsign==0 &&";

  kineCut = OS+accCut+kineCut;

  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooDataSet *datasetMass = (RooDataSet*)fMass->Get("datasetMass");
  RooAddPdf* pdfMASS_Tot = (RooAddPdf*)fMass->Get("pdfMASS_Tot");

  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  ws->import(*datasetMass);
  ws->import(*pdfMASS_Tot);
  RooDataSet *dsAB = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  dsAB->Print();
  cout << "pt: "<<ptLow<<"-"<<ptHigh<<", y: "<<yLow<<"-"<<yHigh<<", Cent: "<<cLow<<"-"<<cHigh<<"%"<<endl;
  cout << "####################################" << endl;
  dsAB->SetName("dsAB");
  ws->import(*dsAB);

  int nBins = min(int( round((ws->var("ctau3DErr")->getMax() - ws->var("ctau3DErr")->getMin())/0.0025) ), 72);
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
  //
  //RooPlot* myPlot_B = ws->var("ctau3DErr"); // modified
  RooPlot* myPlot_B = ws->var("ctau3DErr")->frame(Range(ctauErrLow, ctauErrHigh)); // modified
  TH1D* hTot = (TH1D*)ws->data("dsAB")->createHistogram(("hTot"), *ws->var("ctau3DErr"), 
      Binning(myPlot_B->GetNbinsX(), myPlot_B->GetXaxis()->GetXmin(), myPlot_B->GetXaxis()->GetXmax()));
  //TH1D* hTot = (TH1D*)ws->data("dsAB")->createHistogram(("hTot"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrLow,ctauErrHigh));
  //double ctauErrMin = hTot->GetBinLowEdge(hTot->FindFirstBinAbove(1e-3,1));
  double ctauErrMin;
  //double ctauErrMax = hTot->GetBinCenter(hTot->FindLastBinAbove(1e-3,1));
  double ctauErrMax;
  //cout<<"NBins: "<<hTot->GetNbinsX()<<endl;
  
  cout<<"\n Search Min"<<endl;
  for(int i=0; i<hTot->GetNbinsX()/2; i++){
      cout<<"Bin: "<<i<<", CtauError: "<<hTot->GetBinLowEdge(i)<<", Events: "<<hTot->GetBinContent(i)<<endl;
    if(hTot->GetBinContent(i)<=0&&hTot->GetBinContent(i+1)>=1&&hTot->GetBinContent(i+2)>=1&&hTot->GetBinContent(i+3)>=1){
      cout<<"##### i: "<<i+1<<", CtauError:  "<<hTot->GetBinLowEdge(i+1)<<hTot->GetBinContent(i+1)<<", Events: "<<hTot->GetBinContent(i+1)<<endl;
      cout<<"##### i(Min): "<<i+2<<", CtauError:  "<<hTot->GetBinLowEdge(i+2)<<hTot->GetBinContent(i+2)<<", Events: "<<hTot->GetBinContent(i+2)<<endl;
      cout<<"##### i: "<<i+3<<", CtauError:  "<<hTot->GetBinLowEdge(i+3)<<hTot->GetBinContent(i+3)<<", Events: "<<hTot->GetBinContent(i+3)<<endl;
      ctauErrMin = hTot->GetBinLowEdge(i+1); 
      break;}
  }

  cout<<"\n Search Max"<<endl;
  for(int i=20; i<hTot->GetNbinsX(); i++){
      cout<<"Bin: "<<i<<", CtauError: "<<hTot->GetBinLowEdge(i)+hTot->GetBinWidth(i)<<", Events: "<<hTot->GetBinContent(i)<<endl;
    //if(hTot->GetBinContent(i)>=0&&hTot->GetBinContent(i+1)<=0&&hTot->GetBinContent(i+2)<=0){
    if(hTot->GetBinContent(i)>=0&&hTot->GetBinContent(i+1)<=0){
      cout<<"##### i(Max): "<<i<<", CtauError:  "<<hTot->GetBinLowEdge(i)+hTot->GetBinWidth(i)<<hTot->GetBinContent(i)<<", Events: "<<hTot->GetBinContent(i)<<endl;
      cout<<"##### i: "<<i+1<<", CtauError:  "<<hTot->GetBinLowEdge(i+1)+hTot->GetBinWidth(i+1)<<hTot->GetBinContent(i+1)<<", Events: "<<hTot->GetBinContent(i+1)<<endl;
      cout<<"##### i: "<<i+2<<", CtauError:  "<<hTot->GetBinLowEdge(i+2)+hTot->GetBinWidth(i+2)<<hTot->GetBinContent(i+2)<<", Events: "<<hTot->GetBinContent(i+2)<<endl;
      ctauErrMax = hTot->GetBinLowEdge(i-1)+hTot->GetBinWidth(i-1);
      break;}
  }

  double minRange = (double)(floor(ctauErrMin*100.)/100.);
  double maxRange = (double)(ceil(ctauErrMax*100.)/100.);

  ws->var("ctau3DErr")->setRange("ctauErrWindow", ctauErrMin, ctauErrMax);
  myPlot_B = ws->var("ctau3DErr")->frame(Range(minRange, maxRange)); // modified

  hTot->Draw();
  TH1D* hTot_M = (TH1D*)ws->data("dsAB")->createHistogram(("hTot_M"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrMin,ctauErrMax));
  RooDataHist* totHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hTot_M);
  RooHistPdf* pdfCTAUERR_Tot = new RooHistPdf("pdfCTAUERR_Tot","hist pdf", *ws->var("ctau3DErr"), *totHist);

  //bkg
  RooDataSet* dataw_Bkg = new RooDataSet("dataw_Bkg","TMP_BKG_DATA", (RooDataSet*)ws->data("dataset_SPLOT"),
      RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Bkg_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Bkg_sw");
  TH1D* hBkg = (TH1D*)dataw_Bkg->createHistogram(("hBkg"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrMin,ctauErrMax));
  RooDataHist* bkgHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hBkg);
  RooHistPdf* pdfCTAUERR_Bkg = new RooHistPdf("pdfCTAUERR_Bkg","hist pdf", *ws->var("ctau3DErr"), *bkgHist);
  //data
  RooDataSet* dataw_Sig = new RooDataSet("dataw_Sig","TMP_SIG_DATA", (RooDataSet*)ws->data("dataset_SPLOT"),
      RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Jpsi_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Jpsi_sw");
  TH1D* hSig = (TH1D*)dataw_Sig->createHistogram(("hSig"), *ws->var("ctau3DErr"),Binning(nBins,ctauErrMin,ctauErrMax));
  RooDataHist* sigHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hSig);
  RooHistPdf* pdfCTAUERR_Jpsi = new RooHistPdf("pdfCTAUERR_Jpsi","hist pdf", *ws->var("ctau3DErr"), *sigHist);
  //import
  ws->import(*dataw_Sig);
  ws->import(*dataw_Bkg);
  ws->import(*pdfCTAUERR_Tot);
  ws->import(*pdfCTAUERR_Jpsi);
  ws->import(*pdfCTAUERR_Bkg);

  ws->pdf("pdfCTAUERR_Tot")->setNormRange("ctauErrWindow");//Normalization
  ws->pdf("pdfCTAUERR_Bkg")->setNormRange("ctauErrWindow");
  ws->pdf("pdfCTAUERR_Jpsi")->setNormRange("ctauErrWindow");
  
  cout<<ws->data("dsAB")->numEntries()<<endl;
  cout<<ws->data("dataset_SPLOT")->numEntries()<<endl;

  TCanvas* c_B =  new TCanvas("canvas_B","My plots",554,4,550,520);
  c_B->cd();
  TPad *pad_B_1 = new TPad("pad_B_1", "pad_B_1", 0, 0.16, 0.98, 1.0);
  pad_B_1->SetTicks(1,1);
  pad_B_1->Draw(); pad_B_1->cd();
  //RooPlot* myPlot_B = ws->var("ctau3DErr")->frame(nBins); // bins
  //RooPlot* myPlot_B = ws->var("ctau3DErr")->frame(Range(minRange, maxRange)); // modified
  RooBinning bins(nBins, ctauErrMin, ctauErrMax);
  myPlot_B->SetTitle("");


  c_B->cd();
  c_B->SetLogy();

  pad_B_1->cd();
  gPad->SetLogy();
  RooPlot* myPlot2_B = (RooPlot*)myPlot_B->Clone();
  ws->data("dataset_SPLOT")->plotOn(myPlot2_B,Name("dataCTAUERR_Tot"), MarkerSize(.7),  Binning(bins), DataError(RooAbsData::SumW2));
  ws->pdf("pdfCTAUERR_Tot")->plotOn(myPlot2_B,Name("pdfCTAUERR_Tot"), LineColor(kGreen+1), NormRange("ctauErrWindow"), LineWidth(2));
  ws->data("dataw_Sig")->plotOn(myPlot2_B,Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2), Binning(bins), DataError(RooAbsData::SumW2));
  ws->pdf("pdfCTAUERR_Jpsi")->plotOn(myPlot2_B,Name("pdfCTAUERR_Jpsi"),LineColor(kRed+2), LineWidth(2), NormRange("ctauErrWindow"));
  ws->data("dataw_Bkg")->plotOn(myPlot2_B,Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue+2), MarkerColor(kBlue+2), Binning(bins), DataError(RooAbsData::SumW2));
  ws->pdf("pdfCTAUERR_Bkg")->plotOn(myPlot2_B,Name("pdfCTAUERR_Bkg"), LineColor(kBlue+2), LineWidth(2), NormRange("ctauErrWindow"));
  Double_t YMax = hTot->GetBinContent(hTot->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i=1; i<=hTot->GetNbinsX(); i++) if (hTot->GetBinContent(i)>0) YMin = min(YMin, hTot->GetBinContent(i));
  Double_t Yup(0.),Ydown(0.);
  Yup = YMax*TMath::Power((YMax/YMin), (0.4/(1.0-0.4-0.3)));
  Ydown = YMin/(TMath::Power((YMax/YMin), (0.3/(1.0-0.4-0.3))));
  myPlot2_B->GetYaxis()->SetRangeUser(Ydown,Yup);
  cout<<ws->var("ctau3DErr")->getMin()<<", "<<ws->var("ctau3DErr")->getMax()<<endl;
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
  Double_t outTot = ws->data("dsAB")->numEntries();
  //Double_t outErr = ws->data("dsAB")->reduce(Form("(ctauErr>=%.6f || ctauErr<=%.6f)", ctauErrMax, ctauErrMin))->numEntries();
  Double_t outErr = ws->data("dsAB")->reduce(Form("(ctau3DErr>=%.6f || ctau3DErr<=%.6f)", ctauErrMax, ctauErrMin))->numEntries();
  cout<<"lost evt: ("<<(outErr*100)/outTot<<")%, "<<outErr<<"evts"<<endl;
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
  drawText(Form("Loss: (%.4f%s) %.f evts", outErr*100/outTot, "%", outErr),text_x,text_y-y_diff*3,text_color,text_size);
  //drawText(Form("n_{J/#psi} = %.f #pm %.f",sData.GetYieldFromSWeight("N_Jpsi"),ws->var("N_Jpsi")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
  //drawText(Form("n_{Bkg} = %.f #pm %.f",sData.GetYieldFromSWeight("N_Bkg"), ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
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
  RooHist* hpull_B = frameTMP_B->pullHist("dataCTAUERR_Tot","pdfCTAUERR_Tot",true);
  hpull_B->SetMarkerSize(0.8);
  RooPlot* pullFrame_B = ws->var("ctau3DErr")->frame(Title("Pull Distribution"), Range(minRange, maxRange)) ;
  //RooPlot* pullFrame_B = ws->var("ctau3DErr")->frame(Title("Pull Distribution"), Bins(nCtauResBins), Range(ctauErrLow, ctauErrHigh)) ;
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

  TLine *lB = new TLine(minRange,0, maxRange,0);
  //TLine *lB = new TLine(ctauErrLow, 0, ctauErrHigh, 0);
  lB->SetLineStyle(1);
  lB->Draw("same");
  pad_B_2->Update();
  cout << endl << "************** Finished SPLOT *****************" << endl << endl;

  c_B->Update();
  c_B->SaveAs(Form("../../figs/2DFit_%s/CtauErr/ctauErr_%s_%s.pdf",DATE.Data(),bCont.Data(),kineLabel.Data()));
  
  TFile *outFile = new TFile(Form("../../roots/2DFit_%s/CtauErr/CtauErrResult_%s_%s.root",DATE.Data(),bCont.Data(),kineLabel.Data()),"recreate");
  dataw_Bkg->Write();
  dataw_Sig->Write();
  pdfCTAUERR_Tot->Write();
  pdfCTAUERR_Bkg->Write();
  pdfCTAUERR_Jpsi->Write();
  outFile->Close();

  cout << "Min : " << ctauErrMin << " " << hTot->GetBinLowEdge(ctauErrMin) << endl;
  cout << "Max : " << ctauErrMax << " " << hTot->GetBinLowEdge(ctauErrMax) << endl;
  cout << nBins<< endl;

}
