#include "../../commonUtility.h"
#include "TText.h"
#include "TLine.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TFile.h"
#include "TROOT.h"
#include "../../cutsAndBin.h"
#include <TGraphErrors.h>

using namespace std;

double getSysPR(TString fileName="Dir1", TString fileName2="Dir2");
double getSysNP(TString fileName="Dir1", TString fileName2="Dir2");

void compute_sys_diff_0_180()
{
  const int nPtBinPR = 7;
  double PtBinArrPR[nPtBinPR+1] = {3,4.5,6.5,7.5,9,12,15,50};
  const int nPtBinNP = 5;
  double PtBinArrNP[nPtBinNP+1] = {3,6.5,9,12,15,50};
  //const int nCentBin = 6;
  //double CentBinArr[nCentBin+1] = {0,10,20,30,40,50,90};
  const int nCentBin = 5;
  double CentBinArrPR[nCentBin+1] = {0,20,40,60,80,100};
  double CentBinArrNP[nCentBin+1] = {0,20,40,60,80,100};

  //TH1D* hCentSys_PR[nPtBinPR];
  //TH1D* hCentSys_NP[nPtBinNP];
  //
  //for(int ib=0; ib<nPtBinPR; ib++){
  //  hCentSys_PR[ib] = new TH1D(Form("hCentSys_PR_%.f%.f",PtBinArrPR[ib],PtBinArrPR[ib+1]),";p_{T};Uncertainty",nPtBinPR,PtBinArrPR);
  //}
  //for(int ib=0; ib<nPtBinNP; ib++){
  //  hCentSys_NP[ib] = new TH1D(Form("hCentSys_NP_%.f%.f",PtBinArrNP[ib],PtBinArrNP[ib+1]),";p_{T};Uncertainty",nPtBinNP,PtBinArrNP);
  //}
  TH1D* hCentSys_PR;
  TH1D* hCentSys_NP;
  hCentSys_PR = new TH1D("hCentSys_PR",";p_{T};Uncertainty",nCentBin,CentBinArrPR);
  hCentSys_NP = new TH1D("hCentSys_NP",";p_{T};Uncertainty",nCentBin,CentBinArrNP);
  
  //TString dirNom = "/Users/goni/Documents/CMS/Analysis/Jpsi_Flow/Usercode/Jpsi_v2_PbPb2018/finalPlot/linearfit/";
  TString dirNom = "/Users/goni/Documents/CMS/Analysis/Jpsi_Flow/Usercode/Jpsi_v2_PbPb2018/finalPlot/finalPlot_210929/roots/";
  TString dirAlt = "/Users/goni/Documents/CMS/Analysis/Jpsi_Flow/Usercode/Jpsi_v2_PbPb2018/Systematic/SignalPDFvariation/";

  double v2_sys = 0;

  TString fileNomPR, fileNomNP;
  TString fileAltPR, fileAltNP;
  for(int ib=0; ib<nCentBin; ib++){
    cout<<"PR Bin: "<<ib<<endl;
    fileNomPR = dirNom+Form("v2_pt6.5-50.0_y0.0-2.4_muPt0.0_centrality%.0f-%.0f.root",CentBinArrPR[ib],CentBinArrPR[ib+1]);
    fileAltPR = dirAlt+Form("SignalPDFvariation_v2_pt6.5-50.0_y0.0-2.4_muPt0.0_centrality%.0f-%.0f.root",CentBinArrPR[ib],CentBinArrPR[ib+1]);
    v2_sys = getSysPR(fileNomPR,fileAltPR);
    //hCentSys_PR[ib]->SetBinContent(ib+1, v2_sys); 
    hCentSys_PR->SetBinContent(ib+1, v2_sys); 
    v2_sys=0;
  }

  cout<<""<<endl;
  for(int ib=0; ib<nCentBin; ib++){
    cout<<"NP Bin: "<<ib<<endl;
    fileNomNP = dirNom+Form("v2_pt6.5-50.0_y0.0-2.4_muPt0.0_centrality%.0f-%.0f.root",CentBinArrNP[ib],CentBinArrNP[ib+1]);
    fileAltNP = dirAlt+Form("SignalPDFvariation_v2_pt6.5-50.0_y0.0-2.4_muPt0.0_centrality%.0f-%.0f.root",CentBinArrNP[ib],CentBinArrNP[ib+1]);
    v2_sys = getSysNP(fileNomNP,fileAltNP);
    //hCentSys_NP[ib]->SetBinContent(ib+1, v2_sys); 
    hCentSys_NP->SetBinContent(ib+1, v2_sys); 
    v2_sys=0;
  }

  TFile *wf = new TFile("SignalPDFvariation_sys_0_180.root","recreate");
  wf->cd();
  hCentSys_PR->Write();
  hCentSys_NP->Write();
  //for(int ib=0; ib<nCentBin; ib++){
  //  hCentSys_PR[ib]->Write();
  //}
  //for(int ib=0; ib<nCentBin; ib++){
  //  hCentSys_NP[ib]->Write();
  //}
  //wf->Close();
}

double getSysPR(TString fileName, TString fileName2)
{ 
  double v2_nom = -999;
  double v2_alt = -999;
  double sys = -999; 
  TString fin_nom;
  fin_nom = Form("%s",fileName.Data());
  TFile* f1 = new TFile(fin_nom.Data() );
  
  TH1D *h1 = (TH1D*) f1->Get("prv2");
  cout<<h1->GetBinContent(0)<<endl;
  v2_nom=h1->GetBinContent(0);
  if(v2_nom == -999){cout << "ERROR!!" << endl; return 1;}

  TString fin_alt;
  fin_alt = Form("%s",fileName2.Data());
  TFile* f2 = new TFile(fin_alt.Data() );
  TH1D *h2 = (TH1D*) f2->Get("prv2");
  cout<<h2->GetBinContent(0)<<endl;
  v2_alt=h2->GetBinContent(0);
  //v2_alt=0.05;
  if(v2_nom == -999 || v2_alt == -999){cout << "ERROR!!" << endl; return 1;}
  
  sys = fabs(v2_alt-v2_nom);
  cout << "Unc : " << Form("%.3f",100*double(sys)) << "%%" << endl;
  return sys;
}

double getSysNP(TString fileName, TString fileName2)
{ 
  double v2_nom = -999;
  double v2_alt = -999;
  double sys = -999; 
  TString fin_nom;
  fin_nom = Form("%s",fileName.Data());
  TFile* f1 = new TFile(fin_nom.Data() );
  
  TH1D *h1 = (TH1D*) f1->Get("npv2");
  v2_nom=h1->GetBinContent(0);
  if(v2_nom == -999){cout << "ERROR!!" << endl; return 1;}

  TString fin_alt;
  fin_alt = Form("%s",fileName2.Data());
  TFile* f2 = new TFile(fin_alt.Data() );
  TH1D *h2 = (TH1D*) f2->Get("npv2");
  v2_alt=h2->GetBinContent(0);
  if(v2_nom == -999 || v2_alt == -999){cout << "ERROR!!" << endl; return 1;}
  
  sys = fabs(v2_alt-v2_nom);
  cout << "Unc : " << Form("%.3f",100*double(sys)) << "%%" << endl;
  return sys;
}
