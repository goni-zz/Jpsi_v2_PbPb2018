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

void compute_sys_diff_20_120()
{
  const int nPtBinPR = 7;
  double PtBinArrPR[nPtBinPR+1] = {3,4.5,6.5,7.5,9,12,15,50};
  const int nPtBinNP = 5;
  double PtBinArrNP[nPtBinNP+1] = {3,6.5,9,12,15,50};
  const int nCentBin = 6;
  double CentBinArr[nCentBin+1] = {0,10,20,30,40,50,90};

  //TH1D* hPtSys_PR[nPtBinPR];
  //TH1D* hPtSys_NP[nPtBinNP];
  //
  //for(int ib=0; ib<nPtBinPR; ib++){
  //  hPtSys_PR[ib] = new TH1D(Form("hPtSys_PR_%.f%.f",PtBinArrPR[ib],PtBinArrPR[ib+1]),";p_{T};Uncertainty",nPtBinPR,PtBinArrPR);
  //}
  //for(int ib=0; ib<nPtBinNP; ib++){
  //  hPtSys_NP[ib] = new TH1D(Form("hPtSys_NP_%.f%.f",PtBinArrNP[ib],PtBinArrNP[ib+1]),";p_{T};Uncertainty",nPtBinNP,PtBinArrNP);
  //}
  TH1D* hPtSys_PR;
  TH1D* hPtSys_NP;
  hPtSys_PR = new TH1D("hPtSys_PR",";p_{T};Uncertainty",nPtBinPR,PtBinArrPR);
  hPtSys_NP = new TH1D("hPtSys_NP",";p_{T};Uncertainty",nPtBinNP,PtBinArrNP);
  
  TString dirNom = "/Users/goni/Documents/CMS/Analysis/Jpsi_Flow/Usercode/Jpsi_v2_PbPb2018/finalPlot/finalPlot_210929/roots/";
  TString dirAlt = "/Users/goni/Documents/CMS/Analysis/Jpsi_Flow/Usercode/Jpsi_v2_PbPb2018/Systematic/SignalPDFvariation/";

  double v2_sys = 0;

  TString fileNomPR, fileNomNP;
  TString fileAltPR, fileAltNP;
  for(int ib=0; ib<nPtBinPR; ib++){
    cout<<ib<<endl;
    if(ib<2) fileNomPR = dirNom+Form("v2_pt%.1f-%.1f_y1.6-2.4_muPt0.0_centrality20-120.root",PtBinArrPR[ib],PtBinArrPR[ib+1]);
    else fileNomPR = dirNom+Form("v2_pt%.1f-%.1f_y0.0-2.4_muPt0.0_centrality20-120.root",PtBinArrPR[ib],PtBinArrPR[ib+1]);

    if(ib<2) fileAltPR = dirAlt+Form("SignalPDFvariation_v2_pt%.1f-%.1f_y1.6-2.4_muPt0.0_centrality20-120.root",PtBinArrPR[ib],PtBinArrPR[ib+1]);
    else if(ib==1) fileAltPR = dirAlt+Form("SignalPDFvariation_v2_pt%.1f-%.1f_y1.6-2.4_muPt0.0_centrality20-120.root",PtBinArrPR[ib],PtBinArrPR[ib+1]);
    else fileAltPR = dirAlt+Form("SignalPDFvariation_v2_pt%.1f-%.1f_y0.0-2.4_muPt0.0_centrality20-120.root",PtBinArrPR[ib],PtBinArrPR[ib+1]);
    v2_sys = getSysPR(fileNomPR,fileAltPR);
    //hPtSys_PR[ib]->SetBinContent(ib+1, v2_sys); 
    hPtSys_PR->SetBinContent(ib+1, v2_sys); 
    v2_sys=0;
  }

  for(int ib=0; ib<nPtBinNP; ib++){
    if(ib<1) fileNomNP = dirNom+Form("v2_pt%.1f-%.1f_y1.6-2.4_muPt0.0_centrality20-120.root",PtBinArrNP[ib],PtBinArrNP[ib+1]);
    else fileNomNP = dirNom+Form("v2_pt%.1f-%.1f_y0.0-2.4_muPt0.0_centrality20-120.root",PtBinArrNP[ib],PtBinArrNP[ib+1]);

    if(ib<1) fileAltNP = dirAlt+Form("SignalPDFvariation_v2_pt%.1f-%.1f_y1.6-2.4_muPt0.0_centrality20-120.root",PtBinArrNP[ib],PtBinArrNP[ib+1]);
    else fileAltNP = dirAlt+Form("SignalPDFvariation_v2_pt%.1f-%.1f_y0.0-2.4_muPt0.0_centrality20-120.root",PtBinArrNP[ib],PtBinArrNP[ib+1]);
    v2_sys = getSysNP(fileNomNP,fileAltNP);
    //hPtSys_NP[ib]->SetBinContent(ib+1, v2_sys); 
    hPtSys_NP->SetBinContent(ib+1, v2_sys); 
    v2_sys=0;
  }

  TFile *wf = new TFile("SignalPDFvariation_sys_10_60.root","recreate");
  wf->cd();
  hPtSys_PR->Write();
  hPtSys_NP->Write();
  //for(int ib=0; ib<nPtBinPR; ib++){
  //  hPtSys_PR[ib]->Write();
  //}
  //for(int ib=0; ib<nPtBinNP; ib++){
  //  hPtSys_NP[ib]->Write();
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
