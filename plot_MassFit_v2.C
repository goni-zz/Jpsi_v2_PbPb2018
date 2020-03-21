#include <Riostream.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TNtuple.h>
#include <TString.h>
#include <TH1D.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TInterpreter.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include "cutsAndBin.h"


void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns, 
        const Int_t rows, const Float_t leftOffset=0.,
        const Float_t bottomOffset=0., 
        const Float_t leftMargin=0.2, 
        const Float_t bottomMargin=0.2,
        const Float_t edge=0.05);
void drawDum(float min, float max, double drawXLabel);
void CalEffErr(TGraph *a, double *b);
TGraphAsymmErrors *getEff(TH1F *h1, TH1F *h2);
TGraphAsymmErrors *calcEff(TH1* h1, TH1* h2);
void formatTGraph(TGraph* a, int b, int c, int d);
void formatTCanv(TCanvas* a);
void formatTLeg(TLegend* a);
void formatTH1F(TH1* a, int b, int c, int d);
void SetPad(TH1 *a);
void SetPad(TH2 *a);
void GetJpsiYield(TH1F *a, double *b, int c);
void GetJpsiYieldGauss(TH1F *a, double *b, int c);
void GetJpsiYieldCbEx(TH1F *a, double *b, int c);

TF1 *v2Fit   = new TF1("v2Fit","[1]*(1+2*[0]*TMath::Cos(2.0*x))",-TMath::PiOver2(),TMath::PiOver2());
TF1 *gauss   = new TF1("gauss","0.015*[0]*TMath::Gaus(x,[1],[2],1)",2.6, 3.5);
TF1 *gaussEx = new TF1("gaussEx","0.015*[0]*TMath::Gaus(x,[1],[2],1)+TMath::Exp([3]-x*[4])",2.6, 3.5);
TF1 *Pol2    = new TF1("Pol2","[0]*x*x+[1]*x+[2]",2.6, 3.5);
TF1 *Expo    = new TF1("Expo","TMath::Exp([0]-x*[1])",2.6, 3.5);

double CrystalBallExpo(double *x, double *par) 
{ 
  double N = par[0]; 
 
  double mean = par[1]; 
  double sigma = par[2]; 
  double alpha = par[3]; 
  double n = par[4]; 
 
  double A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2); 
  double B = n/fabs(alpha) - fabs(alpha); 
 
  // factor 0.050 for bin width to get proper integrated yield 
  if ((x[0]-mean)/sigma>-alpha) 
    return 0.009*N*TMath::Gaus(x[0],mean,sigma,1)+TMath::Exp(par[5]-x[0]*par[6]); 
  else 
    return 0.009*N/(sqrt(TMath::TwoPi())*sigma)*A*pow(B-(x[0]-mean)/sigma,-n)+TMath::Exp(par[5]-x[0]*par[6]); 
} 


double CrystalBall(Double_t *x, Double_t *par)
{
  Double_t N = par[0];
  Double_t width = par[1];

  Double_t mean = par[2];
  Double_t sigma = par[3];
  Double_t alpha = par[4];
  Double_t n = par[5];

  Double_t A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
  Double_t B = n/fabs(alpha) - fabs(alpha);

  if ((x[0]-mean)/sigma>-alpha)
    return N*width*TMath::Gaus(x[0],mean,sigma,1);
  else
    return N/(sqrt(TMath::TwoPi())*sigma)*width*A*pow(B-(x[0]-mean)/sigma,-n);
}

Double_t nPol1(Double_t *x, Double_t *par)
{
    Double_t pol = par[0] + par[1]*x[0];
      return pol;
}

Double_t nPol2(Double_t *x, Double_t *par)
{
    Double_t pol = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
      return pol;
}

Double_t nPol3(Double_t *x, Double_t *par)
{
    Double_t pol = par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
      return pol;
}

Double_t nPol4(Double_t *x, Double_t *par)
{
    Double_t pol = par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
      return pol;
}

Double_t CrystalBallPlusPol2(Double_t *x, Double_t *par)
{
    return CrystalBall(x,par) + nPol2(x,&par[4]);
}
/*
  TF1 *CBeX = new TF1("CBeX",CrystalBallEx,8.0,12.0,6);
  CBeX->SetParNames("Yield (#Jpsi)","Mean","Sigma","#alpha","n","a","b");
  CBeX->SetParameter(1, 9.35);
  CBeX->SetParameter(2, 0.055);
  CBeX->SetParameter(3, 1.0);
  CBeX->SetParameter(4, 20.0);
  CBeX->SetLineStyle(1);
  CBeX->SetLineWidth(1.0);
  CBeX->SetLineColor(kViolet);
*/
// fitMTH == 0, fitting
// fitMTH == 1, counting
void plot_MassFit_v2(int cLow = 20, int cHigh = 120,
float ptLow = 3, float ptHigh = 6.5,
float yLow = 1.6, float yHigh = 2.4,
float ctauCut = 0.0495, float SiMuPtCut = 0){
    int fitMTH = 0;
    gROOT->Macro("~/rootlogon.C");

    gStyle->SetOptFit(0);
    //gStyle->SetTitleFont(62,"xyz");
    //gStyle->SetLabelFont(62,"xyz");
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.11);
    gStyle->SetPadLeftMargin(0.11);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetTitleXOffset(1.0);
    gStyle->SetTitleYOffset(1.0);

	TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh) ;

	TFile *f1 = new TFile("/work2/Oniatree/JPsi/skimmed_file/OniaFlowSkim_JpsiTrig_DBAllPD_AddEP_200309.root","READ");
	TH1F *h1 = new TH1F("h1",";mass;",100,2.3,3.8);
	TTree *tree = (TTree *)f1->Get("mmepevt");
	tree->Draw("mass>>h1",Form("fabs(y)>%0.2f && fabs(y)<%0.2f && pt>%0.2f && pt<%0.2f && cBin>%d && cBin<%d && ctau3D>%0.4f",yLow,yHigh,ptLow,ptHigh,cLow,cHigh,ctauCut),"PE");

/*
    TFile *fIn = new TFile("./Hist_Jpsi.root","READ");
    TH1F *bAnaMass = (TH1F*)fIn->Get("mass");
*/


    double yields[10] = {0.0};
    GetJpsiYieldCbEx(h1,yields,fitMTH);
//    GetJpsiYieldCbEx(bAnaMass,yields,fitMTH);
    //GetJpsiYield(bAnaMass,yields,fitMTH);

	TString perc = "%";

    TLatex *lt1 = new TLatex();
    lt1->SetNDC();
    lt1->SetTextSize(0.03);
    TCanvas *c1 = new TCanvas("c1","",900,800);
    int ents[10] = {0};
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    h1->GetYaxis()->SetTitle("Events/(0.009 GeV/c^{2})");
	h1->GetXaxis()->SetRangeUser(2.6,3.5);
	h1->GetXaxis()->SetTitleOffset(1.1);
	h1->GetYaxis()->SetTitleOffset(1.6);
	h1->GetXaxis()->SetTitleSize(0.03);
	h1->GetYaxis()->SetTitleSize(0.03);
	h1->GetXaxis()->SetLabelSize(0.03);
	h1->GetYaxis()->SetLabelSize(0.03);
    h1->Draw("E");
    lt1->SetTextSize(0.03);
    lt1->DrawLatex(0.65,0.89,"CMS Preliminary");
    lt1->DrawLatex(0.65,0.84,"PbPb  #sqrt{s_{NN}} = 5.02 TeV");
    lt1->SetTextSize(0.03);
    //lt1->DrawLatex(0.2,0.60,"Cent. 0 - 20 %");
    //lt1->DrawLatex(0.2,0.52,"6.5 < p_{T}^{J/#psi} < 30 GeV/c");
    lt1->DrawLatex(0.2,0.82,Form("Yields : %0.2f #pm %0.2f",yields[0], yields[1]));
    lt1->DrawLatex(0.2,0.78,Form("Mean : %0.2f #pm %0.4f",yields[4],yields[6]));
    lt1->DrawLatex(0.2,0.74,Form("#sigma : %0.3f MeV/c^{2}",1000*yields[5]));
    lt1->DrawLatex(0.2,0.70,Form("#chi^{2}/ndof : %0.2f / %0.2f",yields[2], yields[3]));
	lt1->DrawLatex(0.2,0.66,Form("%.2f < p_{T}^{#mu#mu} < %.2f GeV/c",ptLow,ptHigh));
	lt1->DrawLatex(0.2,0.62,Form("%.2f < |y^{#mu#mu}| < %.2f",yLow,yHigh));
	lt1->DrawLatex(0.2,0.58,Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()));
	lt1->DrawLatex(0.2,0.54,Form("L_{J/#psi}>%0.4f",ctauCut));
    lt1->SetTextSize(0.04);
    c1->SaveAs(Form("mass_Jpsi_%s.png",kineLabel.Data()));
    c1->SaveAs(Form("mass_Jpsi_%s.pdf",kineLabel.Data()));


}
//(TH1, color, style, pt, eta, rapidity)
void formatTH1F(TH1* a, int b, int c, int d){
    a->SetLineWidth(2);
    a->SetLineStyle(c);
    a->SetMarkerSize(2);
    a->SetLineColor(b);
    a->SetMarkerColor(b);
    a->GetYaxis()->SetTitle("Single #mu Efficiency");
    if(d == 1){	
        a->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV/c)"); 
        //a->GetXaxis()->SetTitle("p_{T} [GeV/c]"); 
    }else if(d == 2){
        a->GetXaxis()->SetTitle("#eta^{#mu}"); 
    }else if(d == 3){
        a->GetXaxis()->SetTitle("rapidity^{#mu}"); 
    }else if(d == 4){
        a->GetXaxis()->SetTitle("Centrality");
    }else if(d == 5){
        a->GetXaxis()->SetTitle("Centrality (%)");
    }
    a->GetXaxis()->CenterTitle(true);
    //a->GetXaxis()->SetLabelSize(0.05);
    //a->GetXaxis()->SetTitleSize(0.05);
    //a->GetXaxis()->SetTitleOffset(1.1);
    //a->GetYaxis()->SetLabelSize(0.05);
    //a->GetYaxis()->SetTitleSize(0.05);
    //a->GetYaxis()->SetTitleOffset(1.2);

}     

void formatTLeg(TLegend* a){

    a->SetFillStyle(0);  
    a->SetFillColor(0); 
    a->SetBorderSize(0);
    //a->SetTextSize(0.03);
    //a->SetTextFont(63);
    //a->SetTextSize(20);
}

void formatTCanv(TCanvas* a){
    a->SetBorderSize(2);
    a->SetFrameFillColor(0);
    a->cd();
    a->SetGrid(1);
    a->SetTickx();
    a->SetTicky();
}

void formatTGraph(TGraph* a, int b, int c, int d)
{

    a->SetMarkerStyle(c);
    a->SetMarkerColor(b);
    a->SetMarkerSize(1.0);
    a->SetLineColor(b);
    a->SetLineWidth(1);
    a->GetXaxis()->SetLabelSize(0.05);
    a->GetXaxis()->SetTitleSize(0.06);
    a->GetXaxis()->SetTitleOffset(1.0);
    a->GetYaxis()->SetTitle("Efficiency");
    //a->GetXaxis()->CenterTitle();
    if(d == 1){	
        a->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
    }else if(d == 2){
        a->GetXaxis()->SetTitle("eta"); 
    }else if(d == 3){
        a->GetXaxis()->SetTitle("rapidity"); 
    }
    a->GetYaxis()->SetRangeUser(0,1);
    a->GetXaxis()->SetRangeUser(0,10);
    a->GetYaxis()->SetLabelSize(0.05);
    a->GetYaxis()->SetTitleSize(0.05);
    a->GetYaxis()->SetTitleOffset(1.3);


}

TGraphAsymmErrors *calcEff(TH1* h1, TH1* h2)
{
    TGraphAsymmErrors *gEfficiency = new TGraphAsymmErrors();
    gEfficiency->BayesDivide(h1,h2);
    return gEfficiency;
}

TGraphAsymmErrors *getEff(TH1F *h1, TH1F *h2)
{

    h1->Sumw2();
    h2->Sumw2();

    TGraphAsymmErrors *result = calcEff(h1,h2);
    return result;
}
void CalEffErr(TGraph *a, double *b){
    const int nbins = 100;
    double x[nbins], y[nbins];
    double sum = 0, errHighSum = 0, errLowSum = 0, sqSumHigh = 0, sqSumLow = 0;
    //double b[3] = 0;
    
    int nBins = a->GetN();
    for(int i=0;i<a->GetN();i++){
        a->GetPoint(i,x[i],y[i]);
        //cout<<"Eff x = "<<x[i]<<" y = "<<y[i]<<endl;
        double eHigh = a->GetErrorYhigh(i);
        double eLow = a->GetErrorYlow(i);
        //cout<<"Err high = "<<eHigh<<" low = "<<eLow<<endl;
        sum += y[i];
        errHighSum += eHigh;
        sqSumHigh += eHigh*eHigh;
        errLowSum += eLow;
        sqSumLow += eLow*eLow;
    }
    b[0] = sum/nBins;
    b[1] = sqrt(sqSumHigh)/nBins;
    b[2] = sqrt(sqSumLow)/nBins;
    //cout<<"b1 : "<<b[0]<<", b2 : "<<b[1]<<", b3 : "<<b[2]<<endl;

    cout<<b[0]<<"^{"<<b[1]<<"}_{"<<b[2]<<"}"<<endl;
    //return b[3];
}

void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge) {
   if (canv==0) {
      Error("makeMultiPanelCanvas","Got null canvas.");
      return;
   }
   canv->Clear();
   
   TPad* pad[columns][rows];

   Float_t Xlow[columns];
   Float_t Xup[columns];
   Float_t Ylow[rows];
   Float_t Yup[rows];
   Float_t PadWidth = 
   (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
   (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
   (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
   (1.0/(1.0-edge))+(Float_t)rows-2.0);
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
      Xlow[i] = Xup[0] + (i-1)*PadWidth;
      Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
      Ylow[i] = Yup[rows-1] + ct*PadHeight;
      Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
      ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
      for(Int_t j=0;j<rows;j++) {
         canv->cd();
         padName = Form("p_%d_%d",i,j);
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
            Xlow[i],Ylow[j],Xup[i],Yup[j]);
         if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
         else pad[i][j]->SetLeftMargin(0);

         if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
         else pad[i][j]->SetRightMargin(0);

         if(j==0) pad[i][j]->SetTopMargin(edge);
         else pad[i][j]->SetTopMargin(0);

         if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
         else pad[i][j]->SetBottomMargin(0);

         pad[i][j]->Draw();
         pad[i][j]->cd();
         pad[i][j]->SetNumber(columns*j+i+1);
      }
   }
}

void drawDum(float min, float max, double drawXLabel){

    TH1D *hdum = new TH1D("hdum","",20,0,1);
    hdum->SetMaximum(max);

    hdum->SetStats(0);

    if(drawXLabel) hdum->SetXTitle("A_{J} = (E_{T}^{j1}-E_{T}^{j2})/(E_{T}^{j1}+E_{T}^{j2})");
    hdum->GetXaxis()->SetLabelSize(20);
    hdum->GetXaxis()->SetLabelFont(43);
    hdum->GetXaxis()->SetTitleSize(22);
    hdum->GetXaxis()->SetTitleFont(43);
    hdum->GetXaxis()->SetTitleOffset(1.5);
    hdum->GetXaxis()->CenterTitle();

    hdum->GetXaxis()->SetNdivisions(905,true);

    hdum->SetYTitle("Ratio");

    hdum->GetYaxis()->SetLabelSize(20);
    hdum->GetYaxis()->SetLabelFont(43);
    hdum->GetYaxis()->SetTitleSize(20);
    hdum->GetYaxis()->SetTitleFont(43);
    hdum->GetYaxis()->SetTitleOffset(2.5);
    hdum->GetYaxis()->CenterTitle();

    hdum->SetAxisRange(0,0.2,"Y");

    hdum->Draw("");

}
void SetPad(TH2 *a){
    a->GetXaxis()->SetLabelSize(16);
    a->GetXaxis()->SetLabelFont(43);
    a->GetXaxis()->SetTitleSize(16);
    a->GetXaxis()->SetTitleFont(43);
    a->GetXaxis()->SetTitleOffset(1.3);
    a->GetXaxis()->CenterTitle();

    a->GetYaxis()->SetLabelSize(16);
    a->GetYaxis()->SetLabelFont(43);
    a->GetYaxis()->SetTitleSize(20);
    a->GetYaxis()->SetTitleFont(43);
    a->GetYaxis()->SetTitleOffset(1.2);
    a->GetYaxis()->CenterTitle();
    //a->SetMaximum(2.0);
    a->SetMinimum(0.0);
}
void SetPad(TH1 *a){
    //a->GetXaxis()->SetLabelSize(20);
    //a->GetXaxis()->SetLabelFont(43);
    //a->GetXaxis()->SetTitleSize(20);
    //a->GetXaxis()->SetTitleSize(16);
    //a->GetXaxis()->SetTitleFont(43);
    //a->GetXaxis()->SetTitleOffset(1.3);
    a->GetXaxis()->CenterTitle();

    //a->GetYaxis()->SetLabelSize(20);
    //a->GetYaxis()->SetLabelFont(43);
    //a->GetYaxis()->SetTitleSize(20);
    //a->GetYaxis()->SetTitleSize(20);
    //a->GetYaxis()->SetTitleFont(43);
    //a->GetYaxis()->SetTitleOffset(1.2);
    a->GetYaxis()->CenterTitle();
    a->SetMaximum(3.0);
    a->SetMinimum(0.0);
}
void GetJpsiYield(TH1F *a, double *b, int c){
    a->GetXaxis()->SetTitle("mass (GeV/c^{2})");
    //a->GetXaxis()->SetTitle("#mu^{-}#mu^{+} mass (GeV/c^{2})");
    //hIncDiMuonMass->Fit("CBeX","LLMR", "", 8.5, 10.3);
    //hIncDiMuonMass->Fit("CBeX","LLMR", "", 8.5, 10.3);
    //

    gauss->SetParameters(1000,3.09,0.010);
    //TF1 *gaussEx = new TF1("gaussEx",CrystalBallExpo,2.6, 3.5);
    TF1 *CbEx = new TF1("CbEx",CrystalBallExpo,2.6, 3.5, 7);
    if(c == 0){
      a->Fit(gauss,"MQ","",2.9, 3.2);
      a->Fit(Expo,"MQ","",2.6, 3.5);
      //gaussEx->SetParameter(0, 9.30);
      gaussEx->SetParameter(0, gauss->GetParameter(0));
      gaussEx->SetParameter(1, gauss->GetParameter(1));
      gaussEx->SetParameter(2, gauss->GetParameter(2));
      gaussEx->SetParameter(3, Expo->GetParameter(0));
      gaussEx->SetParameter(4, Expo->GetParameter(1));
      gaussEx->SetLineStyle(1);
      //gaussEx->SetLineColor(kOrange+7);
      gaussEx->SetLineColor(kViolet+7);
      gaussEx->SetLineWidth(2);
      a->Fit(gaussEx,"LLMQ","",2.6,3.5);
      CbEx->SetParameter(0, gaussEx->GetParameter(0));
      CbEx->SetParameter(1, gaussEx->GetParameter(1));
      CbEx->SetParameter(2, gaussEx->GetParameter(2));
      CbEx->SetParameter(3, 1.5);
      CbEx->SetParameter(4, 2.0);
      CbEx->SetParameter(5, gaussEx->GetParameter(3));
      CbEx->SetParameter(6, gaussEx->GetParameter(4));
      //CbEx->SetLineColor(kRed+2);
      //CbEx->SetLineColor(kOrange+5);
      //CbEx->SetLineColor(kOrange+8);
      CbEx->SetLineColor(kAzure-3);
      //CbEx->SetLineColor(kAzure+7);
      //CbEx->SetLineColor(kOrange-3);
      CbEx->SetLineWidth(3);
      CbEx->SetLineStyle(2);

      a->Fit(CbEx,"LLMQ","",2.6,3.5);
    }
    //a->Fit(gaussEx,"LLMQ","",7.0,11.0);
    //hIncDiMuonMass->Fit("CBeX","LLMR", "", 7.0, 11.0);
    //hIncDiMuonMass->Fit("CBeX","LLMR", "", 9.0, 10.0);
    /*
    b[0] = fabs(CbEx->GetParameter(0));
    b[1] = CbEx->GetParError(0);
    b[2] = CbEx->GetChisquare();
    b[3] = CbEx->GetNDF();
    */
    if(c == 1){
      int Rlow = a->GetXaxis()->FindBin(3.0);
      int Rhigh = a->GetXaxis()->FindBin(3.2);
      b[0] = a->IntegralAndError(Rlow, Rhigh, b[1]);
    }
    if(c == 0){
      b[0] = fabs(CbEx->GetParameter(0));
      cout<<"yields : "<<b[0]<<endl;
      b[1] = CbEx->GetParError(0);
      b[2] = CbEx->GetChisquare();
      b[3] = CbEx->GetNDF();
      b[4] = CbEx->GetParameter(1);
      b[5] = fabs(CbEx->GetParameter(2));
      b[6] = fabs(CbEx->GetParError(1));
    }

}

void GetJpsiYieldGauss(TH1F *a, double *b, int c){
    a->GetXaxis()->SetTitle("mass (GeV/c^{2})");
    //a->GetXaxis()->SetTitle("#mu^{-}#mu^{+} mass (GeV/c^{2})");
    //hIncDiMuonMass->Fit("CBeX","LLMR", "", 8.5, 10.3);
    //hIncDiMuonMass->Fit("CBeX","LLMR", "", 8.5, 10.3);
    //

    gauss->SetParameters(1000,3.09,0.050);
    //TF1 *gaussEx = new TF1("gaussEx",CrystalBallExpo,2.6, 3.5);
    TF1 *CbEx = new TF1("CbEx",CrystalBallExpo,2.6, 3.5, 7);
    if(c == 0){
      a->Fit(gauss,"MQ","",2.9, 3.2);
      a->Fit(Expo,"MQ","",2.6,3.5);
      //gaussEx->SetParameter(0, 9.30);
      gaussEx->SetParameter(0, gauss->GetParameter(0));
      gaussEx->SetParameter(1, gauss->GetParameter(1));
      gaussEx->SetParameter(2, gauss->GetParameter(2));
      gaussEx->SetParameter(3, Expo->GetParameter(0));
      gaussEx->SetParameter(4, Expo->GetParameter(1));
      gaussEx->SetLineStyle(1);
      //gaussEx->SetLineColor(kOrange+7);
      gaussEx->SetLineColor(kViolet+10);
      gaussEx->SetLineWidth(2);
      a->Fit(gaussEx,"LLMQ","",2.6,3.5);
      CbEx->SetParameter(0, gaussEx->GetParameter(0));
      CbEx->SetParameter(1, gaussEx->GetParameter(1));
      CbEx->SetParameter(2, gaussEx->GetParameter(2));
      CbEx->SetParameter(3, 1.5);
      CbEx->SetParameter(4, 2.0);
      CbEx->SetParameter(5, gaussEx->GetParameter(3));
      CbEx->SetParameter(6, gaussEx->GetParameter(4));
      //CbEx->SetLineColor(kRed+2);
      //CbEx->SetLineColor(kOrange+5);
      //CbEx->SetLineColor(kOrange+8);
      CbEx->SetLineColor(kAzure-3);
      //CbEx->SetLineColor(kAzure+7);
      //CbEx->SetLineColor(kOrange-3);
      CbEx->SetLineWidth(3);
      CbEx->SetLineStyle(2);

      //a->Fit(CbEx,"LLMQ","",4.9,6.0);
      a->Fit(gaussEx,"LLMQ","",2.6,3.5);
    }
    //a->Fit(gaussEx,"LLMQ","",7.0,11.0);
    //hIncDiMuonMass->Fit("CBeX","LLMR", "", 7.0, 11.0);
    //hIncDiMuonMass->Fit("CBeX","LLMR", "", 9.0, 10.0);
    /*
    b[0] = fabs(CbEx->GetParameter(0));
    b[1] = CbEx->GetParError(0);
    b[2] = CbEx->GetChisquare();
    b[3] = CbEx->GetNDF();
    */
    if(c == 1){
      int Rlow = a->GetXaxis()->FindBin(3.0);
      int Rhigh = a->GetXaxis()->FindBin(3.2);
      b[0] = a->IntegralAndError(Rlow, Rhigh, b[1]);
    }
    if(c == 0){
      b[0] = fabs(gaussEx->GetParameter(0));
      cout<<"yields : "<<b[0]<<endl;
      b[1] = gaussEx->GetParError(0);
      b[2] = gaussEx->GetChisquare();
      b[3] = gaussEx->GetNDF();
      b[4] = gaussEx->GetParameter(1);
      b[5] = fabs(gaussEx->GetParameter(2));
      b[6] = fabs(gaussEx->GetParError(1));
    }

}


void GetJpsiYieldCbEx(TH1F *a, double *b, int c){
    a->GetXaxis()->SetTitle("mass (GeV/c^{2})");
    //a->GetXaxis()->SetTitle("#mu^{-}#mu^{+} mass (GeV/c^{2})");
    //hIncDiMuonMass->Fit("CBeX","LLMR", "", 8.5, 10.3);
    //hIncDiMuonMass->Fit("CBeX","LLMR", "", 8.5, 10.3);
    //

    gauss->SetParameters(1000,3.09,0.050);
    //TF1 *gaussEx = new TF1("gaussEx",CrystalBallExpo,2.6, 3.5);
    TF1 *CbEx = new TF1("CbEx",CrystalBallExpo,2.6, 3.5, 7);
    if(c == 0){
      a->Fit(gauss,"MQ","",2.9, 3.2);
      a->Fit(Expo,"MQ","",2.6,3.5);
      //gaussEx->SetParameter(0, 9.30);
      gaussEx->SetParameter(0, gauss->GetParameter(0));
      gaussEx->SetParameter(1, gauss->GetParameter(1));
      gaussEx->SetParameter(2, gauss->GetParameter(2));
      gaussEx->SetParameter(3, Expo->GetParameter(0));
      gaussEx->SetParameter(4, Expo->GetParameter(1));
      gaussEx->SetLineStyle(1);
      //gaussEx->SetLineColor(kOrange+7);
      gaussEx->SetLineColor(kViolet+10);
      gaussEx->SetLineWidth(2);
      a->Fit(gaussEx,"LLMQ","",2.6,3.5);
      CbEx->SetParameter(0, gaussEx->GetParameter(0));
      CbEx->SetParameter(1, gaussEx->GetParameter(1));
      CbEx->SetParameter(2, gaussEx->GetParameter(2));
      CbEx->SetParameter(3, 1.5);
      CbEx->SetParameter(4, 2.0);
      CbEx->SetParameter(5, gaussEx->GetParameter(3));
      CbEx->SetParameter(6, gaussEx->GetParameter(4));
      //CbEx->SetLineColor(kRed+2);
      //CbEx->SetLineColor(kOrange+5);
      //CbEx->SetLineColor(kOrange+8);
      CbEx->SetLineColor(kAzure-3);
      //CbEx->SetLineColor(kAzure+7);
      //CbEx->SetLineColor(kOrange-3);
      CbEx->SetLineWidth(3);
      CbEx->SetLineStyle(2);

      a->Fit(CbEx,"LLMQ","",2.6,3.5);
    }
    //a->Fit(gaussEx,"LLMQ","",7.0,11.0);
    //hIncDiMuonMass->Fit("CBeX","LLMR", "", 7.0, 11.0);
    //hIncDiMuonMass->Fit("CBeX","LLMR", "", 9.0, 10.0);
    /*
    b[0] = fabs(CbEx->GetParameter(0));
    b[1] = CbEx->GetParError(0);
    b[2] = CbEx->GetChisquare();
    b[3] = CbEx->GetNDF();
    */
    if(c == 1){
      int Rlow = a->GetXaxis()->FindBin(3.0);
      int Rhigh = a->GetXaxis()->FindBin(3.2);
      b[0] = a->IntegralAndError(Rlow, Rhigh, b[1]);
    }
    if(c == 0){
      b[0] = fabs(gaussEx->GetParameter(0));
      cout<<"yields : "<<b[0]<<endl;
      b[1] = gaussEx->GetParError(0);
      b[2] = gaussEx->GetChisquare();
      b[3] = gaussEx->GetNDF();
      b[4] = gaussEx->GetParameter(1);
      b[5] = fabs(gaussEx->GetParameter(2));
      b[6] = fabs(gaussEx->GetParError(1));
    }

}

