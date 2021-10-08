//Headers{{{
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <Fit/Fitter.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Math/WrappedMultiTF1.h>
#include <Math/ChebyshevPol.h>
#include <HFitInterface.h>
#include "../../commonUtility.h"
#include "../../cutsAndBin.h"
#include "../../HiEvtPlaneList.h"
#include "../../Style.h"
#include "../../tdrstyle.C"
#include "../../CMS_lumi_v2mass.C"

//}}}

using namespace std;

const int nParmM = 11;
const int nParmV = 14;
//Int_t iparmass[nParmM] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Int_t iparmass[nParmM] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
//Int_t iparvn[nParmV] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
Int_t iparvn[nParmV] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};

struct GlobalChi2_width
{
  GlobalChi2_width(ROOT::Math::IMultiGenFunction & f1, ROOT::Math::IMultiGenFunction & f2):
    fChi2_1(&f1), fChi2_2(&f2) {}

  Double_t operator() (const double *par) const
  {
    Double_t p1[nParmM];
    for(Int_t i = 0; i < nParmM; i++) p1[i] = par[iparmass[i]];
    Double_t p2[nParmV];
    for(Int_t i = 0; i < nParmV; i++) p2[i] = par[iparvn[i]];
    return (*fChi2_1)(p1) + (*fChi2_2)(p2);
  }
  const ROOT::Math::IMultiGenFunction * fChi2_1;
  const ROOT::Math::IMultiGenFunction * fChi2_2;
};

//totalYield{{{
Double_t TotalYield(Double_t* x, Double_t* par)
{
  //Double_t N1 = par[0]; //Number of Jpsi yield
  //Double_t Nbkg = par[1]; //Nuber of Bkg
  //Double_t mean = par[2]; //Crystall Ball mean
  //Double_t sigma = par[3]; //Crystall Ball sigma
  //Double_t alpha = par[4]; //crystall ball alpha
  //Double_t n = par[5]; //Crystall ball n
  //Double_t ratio = par[6]; //For fraction of Double Crystall ball
  //Double_t frac = par[7]; //Crystall ball f
  //Double_t Bkgmean = par[8]; //Background fit: exp*erf
  //Double_t Bkgsigma = par[9];
  //Double_t Bkgp0 = par[10];
  //Double_t sigma1_2 = sigma*ratio; //

  Double_t N1 = par[0]; //Number of Jpsi yield
  Double_t Nbkg = par[1]; //Nuber of Bkg
  Double_t mean = par[2]; //Crystall Ball mean
  Double_t sigma = par[3]; //Crystall Ball sigma
  Double_t alpha = par[4]; //crystall ball alpha
  Double_t n = par[5]; //Crystall ball n
  Double_t ratio = par[6]; //For fraction of Double Crystall ball
  Double_t frac = par[7]; //Crystall ball f
  Double_t Bkgp0 = par[8];
  Double_t Bkgp1 = par[9];
  Double_t Bkgp2 = par[10];
  Double_t sigma1_2 = sigma*ratio; //


  //t2 > t1
  Double_t JPsi_t1 = (x[0]-mean)/sigma;
  Double_t JPsi_t2 = (x[0]-mean)/sigma1_2;

  if (alpha < 0)
  {

    cout << "ERROR ::: alpha variable negative!!!! " << endl;
    return -1;
  }
  //Crystall ball fucntion
  Double_t absAlpha = TMath::Abs(alpha);
  Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
  Double_t b = n/absAlpha - absAlpha;


  Double_t JPsi_1 = -1;
  Double_t JPsi_2 = -1;

  if(JPsi_t1 > -alpha){
    JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
  }
  if(JPsi_t1 <= -alpha){
    JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
  }

  if(JPsi_t2 > -alpha){
    JPsi_2 = exp(-JPsi_t2*JPsi_t2/2.);
  }
  if(JPsi_t2 <= -alpha){
    JPsi_2 = a*TMath::Power((b-JPsi_t2),-n);
  }

  double shx = (10*x[0]-37)/3;
  return N1*(frac*JPsi_1 + (1-frac)*JPsi_2)
    //+ Nbkg*(TMath::Exp(-x[0]/Bkgp0)*(TMath::Erf((x[0]-Bkgmean)/(TMath::Sqrt(2)*Bkgsigma))+1)/2.);
    //+ Nbkg*(TMath::Exp(-x[0]/Bkgp0));
    //+ Nbkg*(ROOT::Math::Chebyshev2(x[0],Bkgp0,Bkgp1,Bkgp2));
    //+ Nbkg*(Bkgp0 + Bkgp1*x[0] + Bkgp2*(2.0*x[0]*x[0] - 1.0));
    + Nbkg*(1+Bkgp0*shx + Bkgp1*(2*shx*shx-1) + Bkgp2*(4*shx*shx*shx-3*shx));
}

//totalYieldSig{{{
Double_t TotalYieldSig(Double_t* x, Double_t* par)
{
  Double_t N1 = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];
  Double_t alpha = par[3];
  Double_t n = par[4];
  Double_t ratio = par[5];
  Double_t frac = par[6];
  Double_t sigma1_2 = sigma*ratio;

  //t2 > t1
  Double_t JPsi_t1 = (x[0]-mean)/sigma;
  Double_t JPsi_t2 = (x[0]-mean)/sigma1_2;

  if (alpha < 0)
  {

    cout << "ERROR ::: alpha variable negative!!!! " << endl;
    return -1;
  }

  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
  Double_t b = n/absAlpha - absAlpha;


  Double_t JPsi_1 = -1;
  Double_t JPsi_2 = -1;

  if(JPsi_t1 > -alpha){
    JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
  }
  else if(JPsi_t1 <= -alpha){
    JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
  }

  if(JPsi_t2 > -alpha){
    JPsi_2 = exp(-JPsi_t2*JPsi_t2/2.);
  }
  else if(JPsi_t2 <= -alpha){
    JPsi_2 = a*TMath::Power((b-JPsi_t2),-n);
  }

  return N1*(frac*JPsi_1 + (1-frac)*JPsi_2);
}

Double_t TotalYieldBkg(Double_t* x, Double_t* par)
{
  Double_t Nbkg = par[0]; //Nuber of Bkg
  Double_t cheb0 = par[1];
  Double_t cheb1 = par[2];
  Double_t cheb2 = par[3];

  double shx = (10*x[0]-37)/3;
  Double_t BkgM = Nbkg*(1+cheb0*shx + cheb1*(2*shx*shx-1) + cheb2*(4*shx*shx*shx-3*shx));

  return BkgM;
}

Double_t Totalvnpol1JPsi(Double_t* x, Double_t* par)
{
  Double_t N1 = par[0];
  Double_t Nbkg = par[1];
  Double_t mean = par[2];
  Double_t sigma = par[3];
  Double_t alpha = par[4];
  Double_t n = par[5];
  Double_t ratio = par[6];
  Double_t frac = par[7];
  Double_t Bkgp0 = par[8];
  Double_t Bkgp1 = par[9];
  Double_t Bkgp2 = par[10];
  Double_t c = par[11];
  Double_t c1 = par[12];
  Double_t c2 = par[13];

  Double_t sigma1S_2 = sigma*ratio;

  //t2 > t1
  Double_t JPsi_t1 = (x[0]-mean)/sigma;
  Double_t JPsi_t2 = (x[0]-mean)/sigma1S_2;
  if (alpha < 0)
  {

    cout << "ERROR ::: alpha variable negative!!!! " << endl;
    return -1;
  }

  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
  Double_t b = n/absAlpha - absAlpha;

  Double_t JPsi_1 = -1;
  Double_t JPsi_2 = -1;

  if(JPsi_t1 > -alpha){
    JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
  }
  else if(JPsi_t1 <= -alpha){
    JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
  }

  if(JPsi_t2 > -alpha){
    JPsi_2 = exp(-JPsi_t2*JPsi_t2/2.);
  }
  else if(JPsi_t2 <= -alpha){
    JPsi_2 = a*TMath::Power((b-JPsi_t2),-n);
  }
  Double_t SigM = N1*(frac*JPsi_1 + (1-frac)*JPsi_2);
  Double_t SigM1s = N1*(JPsi_1);
  //Double_t BkgM = Nbkg*(TMath::Exp(-x[0]/Bkgp0)*(TMath::Erf((x[0]-Bkgmean)/(TMath::Sqrt(2)*Bkgsigma))+1)/2.);
  //Double_t BkgM = Nbkg*(TMath::Exp(-x[0]/Bkgp0));
  //Double_t BkgM = Nbkg*(ROOT::Math::Chebyshev2(x[0],Bkgp0,Bkgp1,Bkgp2));
  //Double_t BkgM = Nbkg*(Bkgp0 + Bkgp1*x[0] + Bkgp2*(2.0*x[0]*x[0] - 1.0));
  double shx = (10*x[0]-37)/3;
  Double_t BkgM = Nbkg*(1+Bkgp0*shx + Bkgp1*(2*shx*shx-1) + Bkgp2*(4*shx*shx*shx-3*shx));
  //return c*(SigM1s/(SigM+BkgM)) + (1 - SigM/(SigM+BkgM))*(c3 + c2*x[0] + c1*x[0]*x[0]);
  return c*(SigM/(SigM+BkgM)) + (1 - SigM/(SigM+BkgM))*(c2 + c1*x[0]);

}

//totalvn pol2 bkg Jpsi{{{
Double_t Totalvnpol2JPsi(Double_t* x, Double_t* par)
{
  //Double_t N1 = par[0];
  //Double_t Nbkg = par[1];
  //Double_t mean = par[2];
  //Double_t sigma = par[3];
  //Double_t alpha = par[4];
  //Double_t n = par[5];
  //Double_t ratio = par[6];
  //Double_t frac = par[7];
  //Double_t Bkgmean = par[8];
  //Double_t Bkgsigma = par[9];
  //Double_t Bkgp0 = par[10];
  //Double_t c = par[11];
  //Double_t c1 = par[12];
  //Double_t c2 = par[13];
  //Double_t c3 = par[14];
  Double_t N1 = par[0];
  Double_t Nbkg = par[1];
  Double_t mean = par[2];
  Double_t sigma = par[3];
  Double_t alpha = par[4];
  Double_t n = par[5];
  Double_t ratio = par[6];
  Double_t frac = par[7];
  Double_t Bkgp0 = par[8];
  Double_t Bkgp1 = par[9];
  Double_t Bkgp2 = par[10];
  Double_t c = par[11];
  Double_t c1 = par[12];
  Double_t c2 = par[13];
  Double_t c3 = par[14];

  Double_t sigma1S_2 = sigma*ratio;

  //t2 > t1
  Double_t JPsi_t1 = (x[0]-mean)/sigma;
  Double_t JPsi_t2 = (x[0]-mean)/sigma1S_2;
  if (alpha < 0)
  {

    cout << "ERROR ::: alpha variable negative!!!! " << endl;
    return -1;
  }

  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
  Double_t b = n/absAlpha - absAlpha;

  Double_t JPsi_1 = -1;
  Double_t JPsi_2 = -1;

  if(JPsi_t1 > -alpha){
    JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
  }
  else if(JPsi_t1 <= -alpha){
    JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
  }

  if(JPsi_t2 > -alpha){
    JPsi_2 = exp(-JPsi_t2*JPsi_t2/2.);
  }
  else if(JPsi_t2 <= -alpha){
    JPsi_2 = a*TMath::Power((b-JPsi_t2),-n);
  }
  Double_t SigM = N1*(frac*JPsi_1 + (1-frac)*JPsi_2);
  Double_t SigM1s = N1*(JPsi_1);
  //Double_t BkgM = Nbkg*(TMath::Exp(-x[0]/Bkgp0)*(TMath::Erf((x[0]-Bkgmean)/(TMath::Sqrt(2)*Bkgsigma))+1)/2.);
  //Double_t BkgM = Nbkg*(TMath::Exp(-x[0]/Bkgp0));
  //Double_t BkgM = Nbkg*(ROOT::Math::Chebyshev2(x[0],Bkgp0,Bkgp1,Bkgp2));
  Double_t BkgM = Nbkg*(Bkgp0 + Bkgp1*x[0] + Bkgp2*(2.0*x[0]*x[0] - 1.0));

  //return c*(SigM1s/(SigM+BkgM)) + (1 - SigM/(SigM+BkgM))*(c3 + c2*x[0] + c1*x[0]*x[0]);
  return c*(SigM/(SigM+BkgM)) + (1 - SigM/(SigM+BkgM))*(c3 + c2*x[0] + c1*x[0]*x[0]);

}
//}}}

//totalvn pol3 bkg Jpsi{{{
Double_t Totalvnpol3JPsi(Double_t* x, Double_t* par)
{
  Double_t N1 = par[0];
  Double_t Nbkg = par[3];
  Double_t mean = par[4];
  Double_t sigma = par[5];
  Double_t alpha = par[6];
  Double_t n = par[7];
  Double_t ratio = par[8];
  Double_t frac = par[9];
  Double_t Bkgmean = par[10];
  Double_t Bkgsigma = par[11];
  Double_t Bkgp0 = par[12];
  Double_t c = par[13];
  Double_t c1 = par[14];
  Double_t c2 = par[15];
  Double_t c3 = par[16];
  Double_t c4 = par[17];
  Double_t sigma1S_2 = sigma*ratio;

  //t2 > t1
  Double_t JPsi_t1 = (x[0]-mean)/sigma;
  Double_t JPsi_t2 = (x[0]-mean)/sigma1S_2;
  if (alpha < 0)
  {

    cout << "ERROR ::: alpha variable negative!!!! " << endl;
    return -1;
  }

  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
  Double_t b = n/absAlpha - absAlpha;

  Double_t JPsi_1 = -1;
  Double_t JPsi_2 = -1;

  if(JPsi_t1 > -alpha){
    JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
  }
  else if(JPsi_t1 <= -alpha){
    JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
  }

  if(JPsi_t2 > -alpha){
    JPsi_2 = exp(-JPsi_t2*JPsi_t2/2.);
  }
  else if(JPsi_t2 <= -alpha){
    JPsi_2 = a*TMath::Power((b-JPsi_t2),-n);
  }

  Double_t SigM = N1*(frac*JPsi_1 + (1-frac)*JPsi_2);
  //Double_t BkgM = Nbkg*(TMath::Exp(-x[0]/Bkgp0)*(TMath::Erf((x[0]-Bkgmean)/(TMath::Sqrt(2)*Bkgsigma))+1)/2.);
  Double_t BkgM = Nbkg*(TMath::Exp(-x[0]/Bkgp0)/2.);

  return c*(SigM/(SigM+BkgM)) + (1 - SigM/(SigM+BkgM))*(c4 + c3*x[0] + c2*x[0]*x[0] + c1*x[0]*x[0]*x[0]);

}
//}}}

//Alpha function{{{
Double_t alphaFunct(Double_t* x, Double_t* par)
{
  Double_t N1 = par[0];
  Double_t Nbkg = par[1];
  Double_t mean = par[2];
  Double_t sigma = par[3];
  Double_t alpha = par[4];
  Double_t n = par[5];
  Double_t ratio = par[6];
  Double_t frac = par[7];
  Double_t cheb0 = par[8];
  Double_t cheb1 = par[9];
  Double_t cheb2 = par[10];
  Double_t sigma1S_2 = sigma*ratio;

  //t2 > t1
  Double_t JPsi_t1 = (x[0]-mean)/sigma;
  Double_t JPsi_t2 = (x[0]-mean)/sigma1S_2;
  if (alpha < 0)
  {

    cout << "ERROR ::: alpha variable negative!!!! " << endl;
    return -1;
  }

  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
  Double_t b = n/absAlpha - absAlpha;

  Double_t JPsi_1 = -1;
  Double_t JPsi_2 = -1;

  if(JPsi_t1 > -alpha){
    JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
  }
  else if(JPsi_t1 <= -alpha){
    JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
  }

  if(JPsi_t2 > -alpha){
    JPsi_2 = exp(-JPsi_t2*JPsi_t2/2.);
  }
  else if(JPsi_t2 <= -alpha){
    JPsi_2 = a*TMath::Power((b-JPsi_t2),-n);
  }
  Double_t SigM = N1*(frac*JPsi_1 + (1-frac)*JPsi_2);
  Double_t SigM1s = N1*(JPsi_1);
  double shx = (10*x[0]-37)/3;
  Double_t BkgM = Nbkg*(1+cheb0*shx + cheb1*(2*shx*shx-1) + cheb2*(4*shx*shx*shx-3*shx));

  return (SigM/(SigM+BkgM));

}

Double_t vnPol1BkgAlpha(Double_t* x, Double_t* par)
{
  Double_t N1 = par[0];
  Double_t Nbkg = par[1];
  Double_t mean = par[2];
  Double_t sigma = par[3];
  Double_t alpha = par[4];
  Double_t n = par[5];
  Double_t ratio = par[6];
  Double_t frac = par[7];
  Double_t cheb0 = par[8];
  Double_t cheb1 = par[9];
  Double_t cheb2 = par[10];
  Double_t c1 = par[11];
  Double_t c2 = par[12];
  Double_t c3 = par[13];

  Double_t sigma1_2 = sigma*ratio;

  //t2 > t1
  Double_t JPsi_t1 = (x[0]-mean)/sigma;
  Double_t JPsi_t2 = (x[0]-mean)/sigma1_2;
  if (alpha < 0)
  {
    cout << "ERROR ::: alpha variable negative!!!! " << endl;
    return -1;
  }

  Double_t absAlpha = fabs((Double_t)alpha);
  Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
  Double_t b = n/absAlpha - absAlpha;

  Double_t JPsi_1 = -1;
  Double_t JPsi_2 = -1;

  if(JPsi_t1 > -alpha){
    JPsi_1 = exp(-JPsi_t1*JPsi_t1/2.);
  }
  else if(JPsi_t1 <= -alpha){
    JPsi_1 = a*TMath::Power((b-JPsi_t1),-n);
  }

  if(JPsi_t2 > -alpha){
    JPsi_2 = exp(-JPsi_t2*JPsi_t2/2.);
  }
  else if(JPsi_t2 <= -alpha){
    JPsi_2 = a*TMath::Power((b-JPsi_t2),-n);
  }
  Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
  Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
  Double_t fN_1 = 1./(sigma*(fC+fD));
  Double_t fN_2 = 1./(sigma1_2*(fC+fD));
  Double_t SigM = N1*(fN_1*frac*JPsi_1 + fN_2*(1-frac)*JPsi_2);

  double shx = (10*x[0]-37)/3;
  Double_t BkgM = Nbkg*(1+cheb0*shx + cheb1*(2*shx*shx-1) + cheb2*(4*shx*shx*shx-3*shx));

  return (1 - SigM/(SigM+BkgM))*(c2 + c1*x[0]);

}

//pol2 bkg{{{
Double_t pol1bkg(Double_t* x, Double_t* par)
{
  Double_t c1 = par[0];
  Double_t c2 = par[1];

  return c1*x[0]+c2;
}

Double_t pol2bkg(Double_t* x, Double_t* par)
{
  Double_t c1 = par[0];
  Double_t c2 = par[1];
  Double_t c3 = par[2];

  return c1*x[0]*x[0]+c2*x[0]+c3;
}
//}}}

//pol3 bkg{{{
Double_t pol3bkg(Double_t* x, Double_t* par)
{
  Double_t c1 = par[0];
  Double_t c2 = par[1];
  Double_t c3 = par[2];
  Double_t c4 = par[3];

  return c1*x[0]*x[0]*x[0]+c2*x[0]*x[0]+c3*x[0] + c4;
}
//}}}

void doSimultaneousV2MassFit_weight_pt6p5_7p5_cent40_80_HFUp(
    int ctauCut=-1,
    float ptLow =  6.5, float ptHigh = 7.5,
    float yLow = 0, float yHigh = 2.4,
    int cLow = 40, int cHigh = 80,
    int weight_PR = 0, //PR : 0, NP : 1
    bool fEffW=true, bool fAccW=true, bool isPtW=true, bool isTnP=true,
    float SiMuPtCut = 0, float massLow = 2.6, float massHigh =3.5, bool dimusign=true, 
    int ibkg_vn_sel = fpol1, int PR=2
    )
{

  TString bCont, cutName;

  if(PR==0) bCont="Prompt";
  else if(PR==1) bCont="NonPrompt";
  else if(PR==2) bCont="Inclusive";

  if (ctauCut==-1) cutName="ctauL";
  else if (ctauCut==0) cutName="ctauC";
  else if (ctauCut==1) cutName="ctauR";

  setTDRStyle();
  gStyle->SetOptFit(0000);
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;

  TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh) ;
  TString dimusignString;
  TString kineLabelTot;
  if(dimusign) dimusignString = "OS";
  else if(!dimusign) dimusignString = "SS";

  kineLabelTot = kineLabel + Form("_m%.1f-%.1f",massLow,massHigh) + "_" + dimusignString;

  TF1* fyieldtot;
  TF1* fvntot;
  Double_t v2;
  Double_t v2e;
  Double_t v2_bkg;

  int nParmV_;
  int nParBkg;
  if(ibkg_vn_sel == fpol1) {nParmV_ = 14; nParBkg = 2;} 
  else if(ibkg_vn_sel == fpol2) {nParmV_ = 15; nParBkg = 3;} 
  else if(ibkg_vn_sel == fpol3) {nParmV_ = 18; nParBkg = 4;}
  else{
    cout << "ERROR ::: No Selection for v2 background function!!!!" << endl;
    return;
  }

  TString DATE = "20_40";
  //TString DATE = "210503";
  //TString DATE = "Corr";
  TFile* fMass; TFile* fFinal; TFile* fMC;
  //Get yield distribution{{{
  //TFile* rf = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_v2/upsilonV2/plots/MassV2_190506/Ups_%s.root",kineLabel.Data()),"read");
  fMass = new TFile(Form("../../Macros/2021_09_14/roots/2DFit_%s/Mass/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), "PR", fEffW, fAccW, 1, 1));
  fFinal = new TFile(Form("../../Macros/2021_09_14/roots/2DFit_%s/Final/2DFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), "PR", fEffW, fAccW, 1, 1));
  fMC = new TFile(Form("../../Macros/2021_09_14/roots_MC_Root618_v2Final/Mass/MassFitResult_%s_PRw_Effw%d_Accw%d_PtW%d_TnP%d_test.root", kineLabel.Data(), fEffW, fAccW, isPtW, isTnP));
  RooWorkspace *ws = new RooWorkspace("workspace");
  RooWorkspace *wsmc = new RooWorkspace("workspaceMC");
  RooDataSet *MC_dataset = (RooDataSet*)fMC->Get("datasetMass");
  RooDataSet *datasetMass = (RooDataSet*)fMass->Get("datasetMass");
  RooAddPdf  *pdfMASS_Tot = (RooAddPdf*)fMass->Get("pdfMASS_Tot");
  wsmc->import(*MC_dataset);
  ws->import(*datasetMass);

  TH1D *Fraction1 = (TH1D*)fFinal->Get("Fraction1");
  TH1D *Fraction2 = (TH1D*)fFinal->Get("Fraction2");
  TH1D *Fraction3 = (TH1D*)fFinal->Get("Fraction3");

  double ctauLow =Fraction1->GetBinLowEdge((double)Fraction1->FindFirstBinAbove(1e-3))+Fraction1->GetBinWidth((double)Fraction1->FindFirstBinAbove(1e-3));
  double ctauHigh=Fraction2->GetBinLowEdge((double)Fraction2->FindFirstBinAbove(1e-3))+Fraction2->GetBinWidth((double)Fraction2->FindFirstBinAbove(1e-3));
  double bfrac1=Fraction1->GetBinContent((double)Fraction1->FindFirstBinAbove(1e-3));
  double bfrac2=Fraction2->GetBinContent((double)Fraction2->FindFirstBinAbove(1e-3));
  double bfrac3=Fraction3->GetBinContent((double)Fraction3->FindFirstBinAbove(1e-3));

  TFile *rf;
  if(ctauCut==0) rf = new TFile(Form("../../roots/HFUp_210928/v2mass_hist/%s/Jpsi_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%.5f_%s.root",
        DATE.Data(),kineLabelTot.Data(),fEffW,fAccW,isPtW,isTnP,ctauLow,ctauHigh,cutName.Data()),"read");
  else if(ctauCut==-1) rf = new TFile(Form("../../roots/HFUp_210928/v2mass_hist/%s/Jpsi_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.root",
        DATE.Data(),kineLabelTot.Data(),fEffW,fAccW,isPtW,isTnP,ctauLow,cutName.Data()),"read");
  else if(ctauCut==1) rf = new TFile(Form("../../roots/HFUp_210928/v2mass_hist/%s/Jpsi_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.root",
        DATE.Data(),kineLabelTot.Data(),fEffW,fAccW,isPtW,isTnP,ctauHigh,cutName.Data()),"read");
  else if(ctauCut==2) rf = new TFile(Form("../../roots/HFUp_210928/v2mass_hist/%s/Jpsi_%s_Eff%d_Acc%d_PtW%d_TnP%d_Inc.root",
        DATE.Data(),kineLabelTot.Data(),fEffW,fAccW,isPtW,isTnP),"read");
  //TFile* rf = new TFile("../../Outputs/makeV2Hist_RD/Jpsi_pt6.5-30.0_y0.0-2.4_muPt0.0_centrality20-120_m2.6-3.5_OS.root","read");
  TH1D* h_v2_SplusB = (TH1D*) rf->Get("h_v2_SplusB");  
  TGraphAsymmErrors* g_mass = (TGraphAsymmErrors*) rf->Get("g_mass");  
  TH1D* h_mass = (TH1D*) rf->Get("h_mass");  

  //define function for simultaneous fitting{{{
  TF1* fmass_total = new TF1("fmass_total", TotalYield, massLow, massHigh, nParmM);
  TF1* fvn_simul;
  if(ibkg_vn_sel == fpol1) fvn_simul = new TF1("fvn_simul", Totalvnpol1JPsi, massLow, massHigh, nParmV_);
  else if(ibkg_vn_sel == fpol2) fvn_simul = new TF1("fvn_simul", Totalvnpol2JPsi, massLow, massHigh, nParmV_);
  else if(ibkg_vn_sel == fpol3) fvn_simul = new TF1("fvn_simul", Totalvnpol3JPsi, massLow, massHigh, nParmV_);
  //}}}

  //combine functions{{{
  fmass_total->SetLineColor(2);
  fmass_total->SetLineWidth(1);

  fvn_simul->SetLineColor(2);
  fvn_simul->SetLineWidth(1);

  ROOT::Math::WrappedMultiTF1 wmass(*fmass_total, 1);
  ROOT::Math::WrappedMultiTF1 wvn(*fvn_simul, 1);

  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange massrange;

  massrange.SetRange(massLow,massHigh);
  ROOT::Fit::BinData datamass(opt, massrange);
  ROOT::Fit::FillData(datamass, g_mass);

  ROOT::Fit::DataRange vnrange;
  vnrange.SetRange(massLow,massHigh);
  ROOT::Fit::BinData datavn(opt, vnrange);
  ROOT::Fit::FillData(datavn, h_v2_SplusB);

  ROOT::Fit::Chi2Function mass_chi2(datamass, wmass);
  ROOT::Fit::Chi2Function vn_chi2(datavn, wvn);

  GlobalChi2_width globalChi2(mass_chi2, vn_chi2);

  ROOT::Fit::Fitter fitter;
  //}}}

  Double_t N1_; Double_t Nbkg_;//Par0, 1
  Double_t N1_min; Double_t Nbkg_min;//Par0, 1

  if(ctauCut==-1){
    N1_min = 0;
    N1_ = 100000;
    Nbkg_min = 0;//Par0, 1
    Nbkg_ = 100000;}//Par0, 1
    //N1_min = 0;
    //N1_ = ws->var("N_Jpsi")->getVal()*0.5;
    //Nbkg_min = 0;//Par0, 1
    //Nbkg_ = ws->var("N_Bkg")->getVal()*0.5;}//Par0, 1
  else if(ctauCut==0){
    N1_min = 0;
    N1_ = 230000;
    Nbkg_min = 0;//Par0, 1
    Nbkg_ =330000;}//Par0, 1
    //N1_min = 0;
    //N1_ = ws->var("N_Jpsi")->getVal()*0.3;
    //Nbkg_min = 0;//Par0, 1
    //Nbkg_ = ws->var("N_Bkg")->getVal()*0.3;}//Par0, 1
  else if(ctauCut==1){
    N1_min = 0;
    N1_ = 230000;
    Nbkg_min = 0;//Par0, 1
    Nbkg_ = 200000;}//Par0, 1

  Double_t mean_ = pdgMass.JPsi;//Par2
  Double_t alpha_; Double_t n_;//Par4, 5
  Double_t Bkgp0_;  Double_t Bkgp1_;  Double_t Bkgp2_;//par8, 9, 10
  Double_t c_;  Double_t c1_;  Double_t c2_;//par11, 12, 13
  //Get fitting parameter{{{
  Double_t sigma_ = 0.04;//Par3
  //if(ctauCut==-1){alpha_ = 2.5; n_ = 1.2;}//par4, 5
  //else if(ctauCut==0){alpha_ = 2.1; n_ = 1.2;}
  //else if(ctauCut==1){alpha_ = 2.1; n_ = 1.3;}
  Double_t ratio_ = 0.35;//par6
  Double_t frac_ = 0.45;//par7
  //if(ctauCut==-1){Bkgp0_=0.1;   Bkgp1_=0.1;   Bkgp2_=0.1;}//par8, 9, 10
  if(ctauCut==-1){Bkgp0_=-0.064908;   Bkgp1_=-0.000670835;   Bkgp2_=0.00622086;}
  else if(ctauCut==0){Bkgp0_=-0.064908;   Bkgp1_=-0.000670835;   Bkgp2_=0.00622086;}
  else if(ctauCut==1){Bkgp0_=-0.064908;   Bkgp1_=-0.000670835;   Bkgp2_=0.00622086;}
  
  //if(ctauCut==-1){c_ = 0.05; c1_ = 0.1; c2_ = 0.1;}//par11, 12, 13
  if(ctauCut==-1){c_ = 0.08; c1_ = -0.01428; c2_ = 0.0284097;}
  else if(ctauCut==0){c_ = 0.04; c1_ = -0.01428; c2_ = 0.0284097;}
  else if(ctauCut==1){c_ = 0.07; c1_ = -0.01428; c2_ = 0.0284097;}

  Double_t par0[nParmV];

  gSystem->mkdir(Form ("figs/v2mass_fit_210929/%s",DATE.Data()), kTRUE);
  gSystem->mkdir(Form("roots/v2mass_fit_210929/%s",DATE.Data()), kTRUE);

  par0[0] = ws->var("N_Jpsi")->getVal();
  par0[1] = ws->var("N_Bkg")->getVal();
  par0[2] = ws->var("m_{J/#Psi}")->getVal();
  par0[3] = ws->var("sigma_1_A")->getVal();
  par0[4] = wsmc->var("alpha_1_A")->getVal();
  par0[5] = wsmc->var("n_1_A")->getVal();
  par0[6] = wsmc->var("x_A")->getVal();
  par0[7] = ws->var("f")->getVal();
  par0[8] = ws->var("sl1")->getVal();
  par0[9] = ws->var("sl2")->getVal();
  
  if(ptLow<=6.5) par0[10] = ws->var("sl3")->getVal();
  else par0[10] = Bkgp2_;
  par0[11] = c_;
  par0[12] = c1_;
  par0[13] = c2_;

  Double_t parLimitLow[nParmV]  = {N1_min, Nbkg_min, mean_-0.1,  0.01,
    1.01, 0.4, 0,  0,  
    -30, -30, -30, 
    -0.2, -30, -30};
  //Double_t parLimitHigh[nParmV] = {par0[0]*0.5, par0[1]*0.31, mean_+0.1,  0.9,//Right
  //Double_t parLimitHigh[nParmV] = {par0[0], par0[1], mean_+0.1,  0.9,//Center
  Double_t parLimitHigh[nParmV] = {N1_, Nbkg_, mean_+0.1,  0.9,//Right
    //par0[4]*1.5, par0[5]*1.5, 1, 1,   
    4.1, 3.1, 1, 1,   
    30, 30, 30, 
    0.2, 30, 30};

  fitter.Config().SetParamsSettings(nParmV_, par0);
  for(int ipar = 0; ipar<nParmV_; ipar++){
    //if(ipar<11){
    //  fitter.Config().ParSettings(ipar).SetLimits(par0[ipar],par0[ipar]);
    //}
    //else fitter.Config().ParSettings(ipar).SetLimits(parLimitLow[ipar],parLimitHigh[ipar]);
    
    //fitter.Config().ParSettings(ipar).SetLimits(parLimitLow[ipar],parLimitHigh[ipar]);
    if(ipar==4||ipar==5||ipar==6) fitter.Config().ParSettings(ipar).Fix();
    else fitter.Config().ParSettings(ipar).SetLimits(parLimitLow[ipar],parLimitHigh[ipar]);
  }
  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");
  //}}}

  fitter.FitFCN(nParmV_, globalChi2, 0, datamass.Size()+datavn.Size(), true);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  cout << "****************************************" << endl;
  cout << "	" << "	Lower Limit" << "	" << "Upper Limit" << endl;
  // Print Parameter and Parameter Limits
  for(int i=0; i<14; i++){
    if(i==4)cout << "alpha [" << i << "]    : " << parLimitLow[i] <<"  "<<
      fvn_simul->GetParameter(4)<<"+/-"<< (double)fvn_simul->GetParError(4)<<" " << parLimitHigh[i] << endl;
    if(i==5){cout << "n [" << i << "]        : " << parLimitLow[i] <<"  "<< 
      fvn_simul->GetParameter(5)<<"+/-"<< (double)fvn_simul->GetParError(5)<<"	" << parLimitHigh[i] << endl;}
    if(i==8){cout << "bkgp0 [" << i << "]    : " << parLimitLow[i] <<"  "<< 
      fvn_simul->GetParameter(8)<<"+/-"<< (double)fvn_simul->GetParError(8)<<"	" << parLimitHigh[i] << endl;}
    if(i==9){cout << "bkgp1 [" << i << "]    : " << parLimitLow[i] <<"  "<< 
      fvn_simul->GetParameter(9)<<"+/-"<< (double)fvn_simul->GetParError(9)<<"	" << parLimitHigh[i] << endl;}
    if(i==10){cout <<"bkgp2 [" << i << "]   : " << parLimitLow[i]<<"  "<< 
      fvn_simul->GetParameter(10)<<"+/-"<< (double)fvn_simul->GetParError(10)<<"  " << parLimitHigh[i] << endl;}
    if(i>=11&&i<14){cout << "v2 par [" << i << "]  : "<< 
      parLimitLow[i] << "	" <<"	" << parLimitHigh[i] << endl;}
  }
  cout << " " << endl;

  /*
     if (PR==2) {TString fitRestxt = Form("../roots/v2mass_fit_210929/SimFitResult_%s_%s.txt",bCont.Data(),kineLabel.Data());}
     else if (PR==0 || PR==1) {TString fitRestxt = Form("../roots/v2mass_fit_210929/SimFitResult_%s_%s_ctau_%.2f.txt",bCont.Data(),kineLabel.Data(),ctauCut);}
     ofstream outputFile(fitRestxt.Data());
     result.Print(outputFile);
     outputFile.close();*/

  //Yield fitting result{{{
  fmass_total->SetFitResult(result, iparmass);
  fmass_total->SetRange(massrange().first, massrange().second);
  g_mass->GetListOfFunctions()->Add(fmass_total);

  int nprm_alpha     = 11;
  //TF1* fyield_bkg = new TF1("fyield_bkg", "[0]*( TMath::Exp(-x/[3])*(TMath::Erf((x-[1])/(TMath::Sqrt(2)*[2]))+1)/2. )", massLow, massHigh);
  //  TF1* fyield_bkg = new TF1("fyield_bkg", "[0]*(TMath::Exp(-x/[1]))", massLow, massHigh);
  //TF1* fyield_bkg = new TF1("fyield_bkg", "[0]*([1]+[2]*x+[3]*(2.0*x*x-1.0))", massLow, massHigh);
  TF1* fyield_bkg = new TF1("fyield_bkg", TotalYieldBkg, massLow, massHigh, 4);
  fyield_bkg->FixParameter(0, fmass_total->GetParameter(1));
  fyield_bkg->FixParameter(1, fmass_total->GetParameter(8));
  fyield_bkg->FixParameter(2, fmass_total->GetParameter(9));
  fyield_bkg->FixParameter(3, fmass_total->GetParameter(10));
  //fyield_bkg->FixParameter(0, 2.4015e+07);
  //fyield_bkg->FixParameter(1, 2.8209e-01);
  //fyield_bkg->FixParameter(2, -4.3675e-03);
  //fyield_bkg->FixParameter(3, 1.4486e-02);
  //fyield_bkg->FixParameter(0, fmass_total->GetParameter(1));
  //fyield_bkg->FixParameter(1, fmass_total->GetParameter(8));//Bkgmean
  //fyield_bkg->FixParameter(2, fmass_total->GetParameter(9));//Bkgsigma
  //fyield_bkg->FixParameter(3, fmass_total->GetParameter(10));//Bkgp0

  g_mass->GetListOfFunctions()->Add(fyield_bkg);
  //}}}

  //vn fitting result{{{
  fvn_simul->SetFitResult(result, iparvn);
  fvn_simul->SetRange(vnrange().first, vnrange().second);

  TF1* fyield_sig = new TF1("fyield_sig",TotalYieldSig, massLow, massHigh, 7);
  fyield_sig->FixParameter(0, fvn_simul->GetParameter(0));
  fyield_sig->FixParameter(1, fvn_simul->GetParameter(2));
  fyield_sig->FixParameter(2, fvn_simul->GetParameter(3));
  fyield_sig->FixParameter(3, fvn_simul->GetParameter(4));
  fyield_sig->FixParameter(4, fvn_simul->GetParameter(5));
  fyield_sig->FixParameter(5, fvn_simul->GetParameter(6));
  fyield_sig->FixParameter(6, fvn_simul->GetParameter(7));

  //g_mass->GetListOfFunctions()->Add(fyield_sig);
  //
  TF1* fAlpha = new TF1("fAlpha",alphaFunct,massLow,massHigh,nprm_alpha);
  int prmid_v2val    = 11;
  for(int iparm=0;iparm<nprm_alpha; iparm++){
    fAlpha->FixParameter(iparm,fvn_simul->GetParameter(iparm));
  }
  h_v2_SplusB->GetListOfFunctions()->Add(fAlpha);

  TF1* fvn_bkg;
  TF1* fvn_bkg_alpha;
  if(ibkg_vn_sel == fpol1){
    fvn_bkg = new TF1("fvn_bkg",pol1bkg, massLow, massHigh, nParBkg);
    //fvn_bkg->FixParameter(0, fvn_simul->GetParameter(12));
    //fvn_bkg->FixParameter(1, fvn_simul->GetParameter(13));
    //fvn_bkg->FixParameter(2, fvn_simul->GetParameter(14));
    fvn_bkg->FixParameter(0, fvn_simul->GetParameter(12));
    fvn_bkg->FixParameter(1, fvn_simul->GetParameter(13));
    fvn_bkg->FixParameter(2, fvn_simul->GetParameter(14));
    fvn_bkg_alpha = new TF1("fvn_bkg_alpha",vnPol1BkgAlpha,massLow,massHigh,nParmV_);
    for(int iparm=0;iparm<nParmV_; iparm++){
      if(iparm<=prmid_v2val) fvn_bkg_alpha->FixParameter(iparm,fvn_simul->GetParameter(iparm));
      else if(iparm>=prmid_v2val) fvn_bkg_alpha->FixParameter(iparm,fvn_simul->GetParameter(iparm+1));
    }
  }
  else if(ibkg_vn_sel == fpol2){
    fvn_bkg = new TF1("fvn_bkg",pol2bkg, massLow, massHigh, nParBkg);
    //fvn_bkg->FixParameter(0, fvn_simul->GetParameter(12));
    //fvn_bkg->FixParameter(1, fvn_simul->GetParameter(13));
    //fvn_bkg->FixParameter(2, fvn_simul->GetParameter(14));
    fvn_bkg->FixParameter(0, fvn_simul->GetParameter(12));
    fvn_bkg->FixParameter(1, fvn_simul->GetParameter(13));
    fvn_bkg->FixParameter(2, fvn_simul->GetParameter(14));
  }
  else if(ibkg_vn_sel == fpol3){
    fvn_bkg = new TF1("fvn_bkg",pol3bkg, massLow, massHigh, nParBkg);
    fvn_bkg->FixParameter(0, fvn_simul->GetParameter(12));
    fvn_bkg->FixParameter(1, fvn_simul->GetParameter(13));
    fvn_bkg->FixParameter(2, fvn_simul->GetParameter(14));
    fvn_bkg->FixParameter(3, fvn_simul->GetParameter(15));
  }

  unsigned int nfpxl = 2000;
  //fvn_bkg->SetNpx(nfpxl);
  fvn_simul->SetNpx(nfpxl);
  fAlpha->SetNpx(nfpxl);
  fvn_bkg_alpha->SetNpx(nfpxl);
  //fyield_bkg->SetNpx(nfpxl);
  fmass_total->SetNpx(nfpxl);

  h_v2_SplusB->GetListOfFunctions()->Add(fvn_simul);
  h_v2_SplusB->GetListOfFunctions()->Add(fvn_bkg);

  v2 = fvn_simul->GetParameter(11);
  v2e = fvn_simul->GetParError(11);

  // Drawing 
  fyield_bkg->SetLineColor(kBlue);
  fyield_bkg->SetLineStyle(kDashed);
  fyield_bkg->SetLineWidth(3);
  //  fyield_bkg->SetFillStyle(3344);
  //  fyield_bkg->SetFillColor(kAzure-9);
  //fyield_sig->SetFillStyle(3004);
  //fyield_sig->SetFillColor(kRed);
  fmass_total->SetLineColor(kBlack);
  fmass_total->SetLineWidth(3);
  fvn_simul->SetLineColor(kRed+2);
  fvn_simul->SetLineWidth(3);
  fvn_bkg->SetLineColor(kRed+2);
  fvn_bkg->SetLineWidth(3);
  fvn_bkg->SetLineStyle(kDashed);
  fAlpha->SetLineStyle(6);
  fAlpha->SetLineColor(kGreen+2);
  fAlpha->SetLineWidth(2);
  fvn_bkg_alpha->SetLineStyle(kDashed);
  fvn_bkg_alpha->SetLineColor(kGreen+2);
  fvn_bkg_alpha->SetLineWidth(2);

  SetHistStyle(h_v2_SplusB,0,0);
  SetGraphStyle2(g_mass,0,0);

  g_mass->SetMarkerSize(1);
  //g_mass->SetMinimum(0);
  g_mass->GetYaxis()->SetTitle("Events/(0.025GeV/c^{2})");
  g_mass->GetXaxis()->SetLimits(massLow,massHigh);
  g_mass->GetXaxis()->SetRangeUser(massLow,massHigh);
  g_mass->GetYaxis()->SetTitleOffset(1.7);
  g_mass->GetYaxis()->SetLabelSize(0.055);
  g_mass->GetXaxis()->SetLabelSize(0.055);
  g_mass->GetXaxis()->SetTitleSize(0.07);
  g_mass->GetYaxis()->SetTitleSize(0.07);
  g_mass->GetYaxis()->SetTitleOffset(1.2);
  g_mass->GetXaxis()->SetNdivisions(510);

  h_v2_SplusB->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
  h_v2_SplusB->GetXaxis()->SetRangeUser(massLow,massHigh);
  h_v2_SplusB->GetXaxis()->SetLimits(massLow,massHigh);
  h_v2_SplusB->GetXaxis()->SetLabelSize(0.055);
  h_v2_SplusB->GetXaxis()->SetTitleSize(0.07);
  h_v2_SplusB->GetXaxis()->SetTitleOffset(1.0);
  h_v2_SplusB->GetXaxis()->SetLabelOffset(0.011);
  h_v2_SplusB->GetYaxis()->SetTitle("v_{2}^{S+B}");
  h_v2_SplusB->GetYaxis()->SetRangeUser(-0.05,0.26);
  h_v2_SplusB->GetYaxis()->SetLabelSize(0.055);
  h_v2_SplusB->GetYaxis()->SetTitleSize(0.07);
  h_v2_SplusB->GetYaxis()->SetTitleOffset(1.2);

  double sizeTick = 12;

  float pos_x = 0.03;
  float pos_x_mass = 0.23;
  float pos_y = 0.76;
  float pos_y_diff = 0.071;
  int text_color = 1;
  float text_size = 18;
  TString perc = "%";

  TLegend *leg1 = new TLegend(0.75,0.55,0.95,0.75);
  leg1->AddEntry(g_mass,"Data","lp");
  leg1->AddEntry(fmass_total,"Fit","l");
  //leg1->AddEntry(fyield_sig,"J/#psi Signal","lf");
  leg1->AddEntry(fyield_bkg,"Background","lf");
  leg1->SetLineColor(kWhite);
  //SetLegendStyle(leg1);

  TCanvas* c_mass_v2 = new TCanvas("c_mass_v2","",590,750);
  TPad* pad1 = new TPad("pad1","pad1",0,0.5,1.0,1.0);
  TPad* pad2 = new TPad("pad2","pad2",0,0.0,1.0,0.5);
  c_mass_v2->cd();
  pad1->SetTicks(1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLeftMargin(0.19);
  pad1->SetTopMargin(0.08);
  pad1->cd();
  pad1->SetLogy();
  //double pad1W = pad1->GetWw()*pad1->GetAbsWNDC();
  //double pad1H = pad1->GetWh()*pad1->GetAbsHNDC();
  //double tickScaleX = (pad1->GetUxmax() - pad1->GetUxmin())/(pad1->GetX2()-pad1->GetX1())*pad1H;
  //g_mass->GetXaxis()->SetTickLength(sizeTick/tickScaleX);   
  g_mass->GetXaxis()->SetTitleSize(0);
  g_mass->GetXaxis()->SetLabelSize(0);
  //g_mass->GetXaxis()->SetNdivisions(510);
  g_mass->Draw("AP");
  Double_t YMax = h_mass->GetBinContent(h_mass->GetMaximumBin());
  // Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBinAbove(0.0)) );
  Double_t YMin = 1e99;
  for (int i=1; i<=h_mass->GetNbinsX(); i++) if (h_mass->GetBinContent(i)>0) YMin = min(YMin, h_mass->GetBinContent(i));
  Double_t Yup,Ydown;
  Ydown = YMin/(TMath::Power((YMax/YMin), (0.1/(1.0-0.1-0.4))));
  Yup = YMax*TMath::Power((YMax/YMin), (0.4/(1.0-0.1-0.4)));
  g_mass->GetYaxis()->SetRangeUser(Ydown,Yup);
  cout<<"Range: "<<YMin<<"-"<<YMax<<endl;
  cout<<"Range: "<<Ydown<<"-"<<Yup<<endl;
  
  leg1->Draw("same");
  //fyield_sig->Draw("same");
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_x_mass,pos_y-pos_y_diff*0,text_color,text_size);
  else if(ptLow!=0) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),pos_x_mass,pos_y-pos_y_diff*0,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff*1,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff*1,text_color,text_size);
  //drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  //drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*1,text_color,text_size);
  drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  if(ctauCut==-1){
    drawText(Form("#font[12]{l}_{J/#psi} < %.4f", ctauLow),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
    drawText(Form("NP J/#psi Frac. %.3f", bfrac1),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);}
  else if(ctauCut==0){
    drawText(Form("%.4f < #font[12]{l}_{J/#psi} < %.4f", ctauLow, ctauHigh),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
    drawText(Form("NP J/#psi Frac. %.3f", bfrac3),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);}
  else if(ctauCut==1){drawText(Form("#font[12]{l}_{J/#psi} > %.4f", ctauHigh),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
    drawText(Form("NP J/#psi Frac. %.3f", bfrac2),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);}
  else if(ctauCut==2){drawText(Form("#font[12]{l}_{J/#psi} Inclusive"),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);}
  
  //c_mass_v2->cd();
  //pad2->SetTicks(1,1);
  pad2->SetTicks();
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.17);
  pad2->SetLeftMargin(0.19);
  pad2->cd();
  double pad2W = pad2->GetWw()*pad2->GetAbsWNDC();
  double pad2H = pad2->GetWh()*pad2->GetAbsHNDC();
  /*
     TGraphErrors* g_v2_SplusB = new TGraphErrors(h_v2_SplusB);
     g_v2_SplusB->GetYaxis()->SetRangeUser(-0.17,0.17);
     g_v2_SplusB->GetYaxis()->SetTitle("v_{2}^{S+B}");
     g_v2_SplusB->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
     g_v2_SplusB->GetYaxis()->SetLabelSize(0.055);
     g_v2_SplusB->GetXaxis()->SetLabelSize(0.055);
     g_v2_SplusB->GetXaxis()->SetTitleSize(0.07);
     g_v2_SplusB->GetYaxis()->SetTitleSize(0.07);
     g_v2_SplusB->GetYaxis()->SetTitleOffset(1.2);
     g_v2_SplusB->GetXaxis()->SetTitleOffset(1.0);
     g_v2_SplusB->GetXaxis()->SetLabelOffset(0.011);
     g_v2_SplusB->GetXaxis()->SetLimits(massLow,massHigh);
     g_v2_SplusB->GetXaxis()->SetRangeUser(massLow,massHigh);
     */
  //h_v2_SplusB->GetXaxis()->SetNdivisions(510);
  //tickScaleX = (pad2->GetUxmax() - pad2->GetUxmin())/(pad2->GetX2()-pad2->GetX1())*pad2H;
  //h_v2_SplusB->GetXaxis()->SetTickLength(sizeTick/tickScaleX);
  h_v2_SplusB->Draw("P");
  jumSun(massLow,0,massHigh,0,1,1);
  drawText(Form("v_{2}^{S} = %.3f #pm %.3f",v2,v2e),pos_x_mass+0.45,pos_y+pos_y_diff*2,text_color,text_size+2);
  CMS_lumi_v2mass(pad1,iPeriod,iPos);  
  pad1->Update();
  pad2->Update();
  c_mass_v2->cd();
  pad1->Draw();
  pad2->Draw();
  //c_mass_v2->Update();
  if (ctauCut==0) c_mass_v2->SaveAs(Form("figs/v2mass_fit_210929/%s/v2Mass_%s_ctau_%.5f_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauLow,ctauHigh,cutName.Data()));
  else if(ctauCut==-1)c_mass_v2->SaveAs(Form("figs/v2mass_fit_210929/%s/v2Mass_%s_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauLow,cutName.Data()));
  else if(ctauCut==1)c_mass_v2->SaveAs(Form("figs/v2mass_fit_210929/%s/v2Mass_%s_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauHigh,cutName.Data()));
  else if(ctauCut==2)c_mass_v2->SaveAs(Form("figs/v2mass_fit_210929/%s/v2Mass_%s_Inc.pdf",DATE.Data(),kineLabel.Data()));
  //if (PR==0||PR==1) {c_mass_v2->SaveAs(Form("../figs/v2mass_fit_210929/v2Mass_%s_%s_ctau_%.2f.pdf",bCont.Data(), kineLabel.Data(),ctauCut));}
  /*
     Double_t xmass[200];
     Double_t pullmass[200];

     Float_t Chi2Yield = 0;
     Int_t NdfYield = g_mass->GetN() - fmass_total->GetNumberFreeParameters();

     for(Int_t ibin = 1; ibin <= g_mass->GetN(); ibin++)
     {
     xmass[ibin] = h_mass->GetBinCenter(ibin);
     pullmass[ibin] = (h_mass->GetBinContent(ibin) - fmass_total->Eval(xmass[ibin]))/h_mass->GetBinError(ibin);
     Chi2Yield += pullmass[ibin]*pullmass[ibin];
     }

     TLatex* lt1 = new TLatex();
     lt1->SetNDC();
     lt1->DrawLatex(0.49, 0.5, Form("Chi2/ndf = %.f/%d", Chi2Yield, NdfYield));

     pad2->SetTopMargin(0);
     pad2->SetBottomMargin(0.17);
     pad2->SetLeftMargin(0.19);
     pad2->SetTicks();
     pad2->cd();
     jumSun(massLow,0,massHigh,0,1,1);
     drawText(Form("v_{2}(S) = %.3f #pm %.3f",v2,v2e),pos_x_mass,pos_y,text_color,text_size);

     h_v2_SplusB->Draw("P");

     Double_t xv2[200];
     Double_t pullv2[200];
     Double_t v2y[200];

     Float_t Chi2v2 = 0;
     Int_t Ndfv2 = h_v2_SplusB->GetNbinsX()-fvn_bkg->GetNumberFreeParameters();
     TGraphErrors *gvn = new TGraphErrors(h_v2_SplusB);
     for(Int_t ibin = 0; ibin < gvn->GetN(); ibin++)
     {
     gvn->GetPoint(ibin, xv2[ibin], v2y[ibin]);
     pullv2[ibin] = (v2y[ibin] - fvn_simul->Eval(xv2[ibin]))/gvn->GetErrorY(ibin);
     if(fabs(pullv2[ibin]) < 100)
     {
     Chi2v2 += pullv2[ibin]*pullv2[ibin];
     }
     }

     lt1->DrawLatex(0.5, 0.38, Form("Chi2/ndf = %.f/%d", Chi2v2, Ndfv2));
  //}}}
  */
  /*
  //Pull Distribution 
  TCanvas* c_pull = new TCanvas("c_pull","",590,750);
  TPad* pad_pull1 = new TPad("pad_pull1","pad_pull1",0,0.5,1.0,1.0);
  TPad* pad_pull2 = new TPad("pad_pull2","pad_pull2",0,0.0,1.0,0.5);
  c_pull->cd();
  pad_pull1->SetTicks(1,1);
  pad_pull1->SetBottomMargin(0);
  pad_pull1->SetLeftMargin(0.19);
  pad_pull1->SetTopMargin(0.08);
  pad_pull1->Draw();
  pad_pull1->cd();

  TH1D* h_mass_pull = (TH1D*) h_mass->Clone("h_mass_pull"); h_mass->Reset();
  for(Int_t ibin = 1; ibin <= h_mass->GetNbinsX(); ibin++)
  {
  xmass[ibin] = h_mass->GetBinCenter(ibin);
  pullmass[ibin] = (h_mass->GetBinContent(ibin) - fmass_total->Eval(xmass[ibin]))/h_mass->GetBinError(ibin);
  Chi2Yield += pullmass[ibin]*pullmass[ibin];
  }

  TLatex* lt1 = new TLatex();
  lt1->SetNDC();
  lt1->DrawLatex(0.49, 0.5, Form("Chi2/ndf = %.f/%d", Chi2Yield, NdfYield));


  Double_t xv2[200];
  Double_t pullv2[200];
  Double_t v2y[200];

  Float_t Chi2v2 = 0;
  Int_t Ndfv2 = h_v2_SplusB->GetNbinsX()-fvn_bkg->GetNumberFreeParameters();
  TGraphErrors *gvn = new TGraphErrors(h_v2_SplusB);
  for(Int_t ibin = 0; ibin < gvn->GetN(); ibin++)
  {
  gvn->GetPoint(ibin, xv2[ibin], v2y[ibin]);
  pullv2[ibin] = (v2y[ibin] - fvn_simul->Eval(xv2[ibin]))/gvn->GetErrorY(ibin);
  if(fabs(pullv2[ibin]) < 100)
  {
  Chi2v2 += pullv2[ibin]*pullv2[ibin];
  }
  }

  lt1->DrawLatex(0.5, 0.38, Form("Chi2/ndf = %.f/%d", Chi2v2, Ndfv2));
  //}}}
  */
  TFile *wf;
  if (ctauCut==0) {wf = new TFile(Form("roots/v2mass_fit_210929/%s/SimFitResult_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%.5f_%s.root", DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLow,ctauHigh,cutName.Data()),"recreate"); wf->cd();}
  else if(ctauCut==-1) {wf = new TFile(Form("roots/v2mass_fit_210929/%s/SimFitResult_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.root", DATE.Data(), kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLow,cutName.Data()),"recreate"); wf->cd();}
  else if(ctauCut==1) {wf = new TFile(Form("roots/v2mass_fit_210929/%s/SimFitResult_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.root", DATE.Data(), kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauHigh,cutName.Data()),"recreate"); wf->cd();}
  //  else if (PR==0 || PR==1) {TFile *wf = new TFile(Form("../roots/v2mass_fit_210929/SimFitResult_%s_%s_ctau_%.2f.root",bCont.Data(), kineLabel.Data(), ctauCut),"recreate"); wf->cd();}
  //store individual function{{{
  fyieldtot = (TF1*) fmass_total->Clone();
  fyieldtot->SetName("massfit");
  fyieldtot->Write();

  fvntot = (TF1*) fvn_simul->Clone();
  fvntot->SetName("vnfit");
  fvntot->Write();

  fyield_bkg->Write();
  fyield_sig->Write();
  fvn_bkg->Write();
  fAlpha->Write();
  fvn_bkg_alpha->Write();
  //}}}

  //v2_bkg[ipt] = fvn_simul->GetParameter(14) + fvn_simul->GetParameter(14) * JPsi_mass;
  //get Chi2{{{
  Double_t ptar = (ptHigh+ptLow)/2;
  TGraphErrors* v2plot = new TGraphErrors();
  v2plot->SetPoint(0,ptar,v2);
  v2plot->SetPointError(0,0,v2e);
  v2plot->SetTitle("");
  v2plot->SetName("v2vspt");
  v2plot->Write();
  }
