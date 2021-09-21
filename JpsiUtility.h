#include <iostream>
#include "RooPlot.h"
#include "rootFitHeaders.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"

using namespace RooFit;
using namespace std;


bool bkgfit=true;
bool sigfit=true;
bool totfit=true;
float massLow = 2.6, massHigh = 3.5;
int   nMassBin  = 36; //(massHigh-massLow)*30;

float ctauLow = -4, ctauHigh = 7.;
float ctauResLow = -10, ctauResHigh = 10;
int   nCtauErrBins = 144;
int   nCtauResBins = 72;
int   nCtauBins  = (ctauHigh-ctauLow)*10;
int   nCtauTrueBins  = 50;

double ctauErrLow = 1e-6, ctauErrHigh = 0.22;
double binWidth = 0.0025;

float text_x = 0.15;
float text_y = 0.816;
float y_diff = 0.05;
float text_size = 13;
int text_color = 1;

int nCPU = 24; // Change to appropriate your cpu

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
    //myws->import(chi2Var); myws->import(ndofVar);
  } else {
    RooRealVar chi2Var((string("chi2_")+pdfLabel).c_str(),(string("chi2_")+pdfLabel).c_str(),chi2);
    RooRealVar ndofVar((string("ndof_")+pdfLabel).c_str(),(string("ndof_")+pdfLabel).c_str(),ndof);
    //myws->import(chi2Var); myws->import(ndofVar);
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
  pullFrame->GetYaxis()->SetRangeUser(-7,7);
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

void getRange(TH1* hist, int nMaxBins, vector<double>& rangeErr)
{
  // 1) Find the bin with the maximum Y value
  int binMaximum = hist->GetMaximumBin();
  // 2) Loop backward and find the first bin
  int binWithContent = -1;
  int firstBin = 1;
  for (int i=binMaximum; i>0; i--) {
    if (hist->GetBinContent(i)>0.0) {
      if ( binWithContent>0 && ((binWithContent-i) > nMaxBins) && hist->GetBinContent(i)<5.0 ) { firstBin = binWithContent; break; }
      else { binWithContent = i; }
    }
  }
  // 3) Loop forward and find the last bin
  binWithContent = -1;
  int lastBin = hist->GetNbinsX();
  for (int i=binMaximum; i<hist->GetNbinsX(); i++) {
    if (hist->GetBinContent(i)>0.0) {
      if ( binWithContent>0 && ((i - binWithContent) > nMaxBins) && hist->GetBinContent(i)<5.0 ) { lastBin = binWithContent+1; break; }
      else { binWithContent = i; }
    }
  }
  // 4) Build the set of bins
  int startBin = ( (firstBin>1) ? (firstBin-1) : firstBin );
  const int nNewBins = lastBin - startBin + 1;
  double binning[nNewBins+1];
  binning[0] = hist->GetXaxis()->GetXmin();
  binning[nNewBins] = hist->GetXaxis()->GetXmax();
  for (int i=1; i<nNewBins; i++) {
    int iBin = startBin + i;
    binning[i] = hist->GetBinLowEdge(iBin);
  }
  rangeErr.push_back(binning[(firstBin>1)?1:0]);
  rangeErr.push_back(binning[nNewBins-1]);

  return;
};
