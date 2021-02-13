#include <iostream>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include <TMath.h>
#include "commonUtility.h"
#include "cutsAndBin.h"
#include "HiEvtPlaneList.h"
#include "Style.h"
#include "tdrstyle.C"
#include "CMS_lumi_v2mass.C"
#include "rootFitHeaders.h"
using namespace std;
using namespace RooFit;

using namespace hi;

double getAccWeight(TH1D* h = 0, double pt = 0);
double getEffWeight(TH1D* h = 0, double pt = 0);
void GetHistSqrt(TH1D* h1 =0, TH1D* h2=0);
void GetHistBkg(TH1D* h1 =0, TH1D* h2=0);

void v2mass_hist_weight(int cLow = 0, int cHigh = 200,
    float ptLow =  12, float ptHigh =  50,
    float yLow = 1.6, float yHigh = 2.4,
    float SiMuPtCut = 0, float massLow = 2.6, float massHigh =3.5, bool dimusign=true, bool fAccW = false, bool fEffW = false, bool isMC = false, int dtype = 1, 
    float ctauCut =0.2,
	int weight_PR = 0, // PR : 0, NP : 1
    int PR=2) //PR 0: PR, 1: NP, 2: Inc.
{
  //Basic Setting
  //gStyle->SetOptStat(0);
  gStyle->SetOptStat(000000000);
  gROOT->ForceStyle();
  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  TH1::SetDefaultSumw2();
  //TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh, 1) ;
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh) ;
  TString dimusignString;
  if(dimusign) dimusignString = "OS";
  else if(!dimusign) dimusignString = "SS";
  TString sample;
  TString bCont;
  TString wName;

  if(PR==0) bCont="Prompt";
  else if(PR==1) bCont="NonPrompt";
  else if(PR==2) bCont="Inclusive";

  if(weight_PR==0) wName="prompt";
  else if(weight_PR==1) wName="nprompt";

  cout<<"Data type: "<<dtype<<endl;

  TChain *tree = new TChain("mmepevt");
  if(!isMC){
    TString f1 = "skimmedFiles/OniaFlowSkim_JpsiTrig_DBAllPD_isMC0_HFNom_201127.root";
    //TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/OniaFlowSkim_JpsiTrig_DBPD_isMC0_HFNom_AddEP_200217.root";
    //TString f2 = "/Users/goni/Downloads/ONIATREESKIMFILE/OniaFlowSkim_JpsiTrig_DBPeriPD_isMC0_HFNom_AddEP_Peri_200217.root";
    tree->Add(f1.Data());
    //tree->Add(f2.Data());
    sample="RD";}
  else if(dtype==1){
    TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/";
    tree->Add(f1.Data());
    sample="MC_PR";}
  else if(dtype==2){
    TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/";
    tree->Add(f1.Data());
    sample="MC_NP";}

  //READ Input Skimmed File
  //TString f1 = "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/skimmedFiles/OniaFlowSkim_JpsiTrig_DBPD_isMC0_HFNom_190716.root";
  //TString f2 = "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/skimmedFiles/OniaFlowSkim_JpsiTrig_DBPsiPeriPD_isMC0_HFNom_190716.root";
  //TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/OniaFlowSkim_JpsiTrig_AllPD_HFNom.root";
  //TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/OniaFlowSkim_MC_Jpsi_HFNom.root";
  //TString f1 = "/Users/goni/Downloads/ONIATREESKIMFILE/OniaFlowSkim_BJpsiTrig_DBPD_isMC1_HFNom_200217.root";

  //Get Correction histograms
  bool isTnP = false;
  bool isPtW = false;
  //  TFile *fEff = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/Efficiency/mc_eff_vs_pt_TnP%d_PtW1_OfficialMC_Y%dS_muPtCut3.5.root",isTnP,state),"read");
  TFile *fEff = new TFile(Form("mc_eff_vs_pt_cent_rap_%s_pbpb_Jpsi.root",wName.Data()),"read");
  TH1D* hEffPt[12];
  hEffPt[0] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to20_absy0_0p6",isTnP,isPtW));
  hEffPt[1] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to20_absy0p6_1p2",isTnP,isPtW));
  hEffPt[2] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to20_absy1p2_1p6",isTnP,isPtW));
  hEffPt[3] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to20_absy1p6_2p4",isTnP,isPtW));
  hEffPt[4] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent20to100_absy0_0p6",isTnP,isPtW));
  hEffPt[5] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent20to100_absy0p6_1p2",isTnP,isPtW));
  hEffPt[6] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent20to100_absy1p2_1p6",isTnP,isPtW));
  hEffPt[7] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent20to100_absy1p6_2p4",isTnP,isPtW));
  hEffPt[8] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent100to180_absy0_0p6",isTnP,isPtW));
  hEffPt[9] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent100to180_absy0p6_1p2",isTnP,isPtW));
  hEffPt[10] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent100to180_absy1p2_1p6",isTnP,isPtW));
  hEffPt[11] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent100to180_absy1p6_2p4",isTnP,isPtW));
  //  TFile *fAcc = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/Acceptance/acceptance_wgt_%dS_pt0_50_20190813_dNdptWeighted.root",state),"read");
  TFile *fAcc = new TFile(Form("mc_Acc_vs_pt_cent_rap_%s_pbpb_Jpsi_PtW1.root",wName.Data()),"read");
  TH1D* hAccPt[12];
  hAccPt[0] = (TH1D*) fAcc -> Get(Form("mc_Acc_vs_pt_PtW%d_cent0to20_absy0_0p6",isPtW));
  hAccPt[1] = (TH1D*) fAcc -> Get(Form("mc_Acc_vs_pt_PtW%d_cent0to20_absy0p6_1p2",isPtW));
  hAccPt[2] = (TH1D*) fAcc -> Get(Form("mc_Acc_vs_pt_PtW%d_cent0to20_absy1p2_1p6",isPtW));
  hAccPt[3] = (TH1D*) fAcc -> Get(Form("mc_Acc_vs_pt_PtW%d_cent0to20_absy1p6_2p4",isPtW));
  hAccPt[4] = (TH1D*) fAcc -> Get(Form("mc_Acc_vs_pt_PtW%d_cent20to100_absy0_0p6",isPtW));
  hAccPt[5] = (TH1D*) fAcc -> Get(Form("mc_Acc_vs_pt_PtW%d_cent20to100_absy0p6_1p2",isPtW));
  hAccPt[6] = (TH1D*) fAcc -> Get(Form("mc_Acc_vs_pt_PtW%d_cent20to100_absy1p2_1p6",isPtW));
  hAccPt[7] = (TH1D*) fAcc -> Get(Form("mc_Acc_vs_pt_PtW%d_cent20to100_absy1p6_2p4",isPtW));
  hAccPt[8] = (TH1D*) fAcc -> Get(Form("mc_Acc_vs_pt_PtW%d_cent100to180_absy0_0p6",isPtW));
  hAccPt[9] = (TH1D*) fAcc -> Get(Form("mc_Acc_vs_pt_PtW%d_cent100to180_absy0p6_1p2",isPtW));
  hAccPt[10] = (TH1D*) fAcc -> Get(Form("mc_Acc_vs_pt_PtW%d_cent100to180_absy1p2_1p6",isPtW));
  hAccPt[11] = (TH1D*) fAcc -> Get(Form("mc_Acc_vs_pt_PtW%d_cent100to180_absy1p6_2p4",isPtW));

  //SetBranchAddress
  const int nMaxDimu = 1000;
  float mass[nMaxDimu];
  float qxa[nMaxDimu];
  float qya[nMaxDimu];
  float qxb[nMaxDimu]; 
  float qyb[nMaxDimu];
  float qxc[nMaxDimu];
  float qyc[nMaxDimu];
  float qxdimu[nMaxDimu];
  float qydimu[nMaxDimu];
//  float P[nMaxDimu];
//  float Px[nMaxDimu];
//  float Py[nMaxDimu];
//  float Pz[nMaxDimu];
  float pt[nMaxDimu];
  float y[nMaxDimu]; 
  float pt1[nMaxDimu];
  float pt2[nMaxDimu];
  float eta[nMaxDimu];
  float eta1[nMaxDimu];
  float eta2[nMaxDimu];
  float ctau3D[nMaxDimu];
  Int_t cBin;
  Int_t event; 
  Int_t nDimu; 
  float vz;
  int recoQQsign[nMaxDimu];
  double weight;

  TBranch *b_event;
  TBranch *b_cBin;
  TBranch *b_nDimu;
  TBranch *b_vz;
  TBranch *b_mass;
  TBranch *b_recoQQsign;
//  TBranch *b_P;
//  TBranch *b_Px;
//  TBranch *b_Py;
//  TBranch *b_Pz;
  TBranch *b_pt;
  TBranch *b_y;
  TBranch *b_eta;
  TBranch *b_eta1;
  TBranch *b_eta2;
  TBranch *b_ctau3D;
  TBranch *b_pt1;
  TBranch *b_pt2;
  TBranch *b_qxa;
  TBranch *b_qxb;
  TBranch *b_qxc;
  TBranch *b_qxdimu;
  TBranch *b_qya;
  TBranch *b_qyb;
  TBranch *b_qyc;
  TBranch *b_qydimu;
  TBranch *b_weight;


  tree -> SetBranchAddress("event", &event, &b_event);
  tree -> SetBranchAddress("cBin", &cBin, &b_cBin);
  tree -> SetBranchAddress("nDimu", &nDimu, &b_nDimu);
  tree -> SetBranchAddress("vz", &vz, &b_vz);
  tree -> SetBranchAddress("recoQQsign", recoQQsign, &b_recoQQsign);
  tree -> SetBranchAddress("mass", mass, &b_mass);
  tree -> SetBranchAddress("y", y, &b_y);
//  tree -> SetBranchAddress("P", P, &b_P);
//  tree -> SetBranchAddress("Px", Px, &b_Px);
//  tree -> SetBranchAddress("Py", Py, &b_Py);
//  tree -> SetBranchAddress("Pz", Pz, &b_Pz);
  tree -> SetBranchAddress("pt", pt, &b_pt);
  tree -> SetBranchAddress("pt1", pt1, &b_pt1);
  tree -> SetBranchAddress("pt2", pt2, &b_pt2);
  tree -> SetBranchAddress("eta", eta, &b_eta);
  tree -> SetBranchAddress("eta1", eta1, &b_eta1);
  tree -> SetBranchAddress("eta2", eta2, &b_eta2);
  tree -> SetBranchAddress("ctau3D", ctau3D, &b_ctau3D);
  tree -> SetBranchAddress("qxa", qxa, &b_qxa);
  tree -> SetBranchAddress("qxb", qxb, &b_qxb);
  tree -> SetBranchAddress("qxc", qxc, &b_qxc);
  tree -> SetBranchAddress("qxdimu", qxdimu, &b_qxdimu);
  tree -> SetBranchAddress("qya", qya, &b_qya);
  tree -> SetBranchAddress("qyb", qyb, &b_qyb);
  tree -> SetBranchAddress("qyc", qyc, &b_qyc);
  tree -> SetBranchAddress("qydimu", qydimu, &b_qydimu);
  tree -> SetBranchAddress("weight", &weight, &b_weight);

  //PR
  //if(PR==0){
  //  if(ptLow==3&&ptHigh==4)        ctauCut=0.0625;
  //  if(ptLow==3&&ptHigh==4.5)      ctauCut=0.0595;
  //  else if(ptLow==4&&ptHigh==5)   ctauCut=0.0525;
  //  else if(ptLow==4.5&&ptHigh==6.5)   ctauCut=0.0485;
  //  else if(ptLow==5&&ptHigh==6.5) ctauCut=0.0455;
  //  else if(ptLow==6.5&&ptHigh==8) ctauCut=0.0405;
  //  else if(ptLow==6.5&&ptHigh==9) ctauCut=0.0395;
  //  else if(ptLow==8&&ptHigh==12)  ctauCut=0.0335;
  //  else if(ptLow==8&&ptHigh==15)  ctauCut=0.0315;
  //  else if(ptLow==9&&ptHigh==12)  ctauCut=0.0325;
  //  else if(ptLow==9&&ptHigh==15)  ctauCut=0.0305;
  //  else if(ptLow==12&&ptHigh==20) ctauCut=0.0265;
  //  else if(ptLow==12&&ptHigh==30) ctauCut=0.0255;
  //  else if(ptLow==15&&ptHigh==30) ctauCut=0.0235;
  //  else if(ptLow==0&&ptHigh==30)  ctauCut=0.0765;
  //  else if(cLow==0&&cHigh==20)    ctauCut=0.0385;
  //  else if(cLow==20&&cHigh==40)   ctauCut=0.0375;
  //  else if(cLow==40&&cHigh==60)   ctauCut=0.0375;
  //  else if(cLow==60&&cHigh==80)   ctauCut=0.0375;
  //  else if(cLow==60&&cHigh==120)  ctauCut=0.0375;
  //  else if(cLow==80&&cHigh==100)  ctauCut=0.0365;
  //  else if(cLow==100&&cHigh==180) ctauCut=0.0365;
  //  else if(cLow==120&&cHigh==180) ctauCut=0.0365;}//All
  ////NP
  //else if(PR==1){
  //  if(ptLow==3.0&&ptHigh==6.5)    ctauCut=0.0975;
  //  else if(ptLow==3&&ptHigh==4)   ctauCut=0.1235;
  //  else if(ptLow==3&&ptHigh==4.5) ctauCut=0.1155;
  //  else if(ptLow==4&&ptHigh==6.5) ctauCut=0.0935;
  //  else if(ptLow==4.5&&ptHigh==6) ctauCut=0.0965;
  //  else if(ptLow==4.5&&ptHigh==6.5) ctauCut=0.0935;
  //  else if(ptLow==6.5&&ptHigh==8) ctauCut=0.0785;
  //  else if(ptLow==6.5&&ptHigh==9) ctauCut=0.0755;
  //  else if(ptLow==8&&ptHigh==12)  ctauCut=0.0655;
  //  else if(ptLow==9&&ptHigh==12)  ctauCut=0.0635;
  //  else if(ptLow==9&&ptHigh==15)  ctauCut=0.0605;
  //  else if(ptLow==12&&ptHigh==30) ctauCut=0.0505;
  //  else if(ptLow==10&&ptHigh==15) ctauCut=0.0585;
  //  else if(ptLow==15&&ptHigh==30) ctauCut=0.0455;
  //  else if(ptLow==15&&ptHigh==50) ctauCut=0.0455;
  //  else if(ptLow==0&&ptHigh==30)  ctauCut=0.0765;}//All
  //crosscheck+tnp+genweight
  //  else if(PR==1){
  //  else if(ptLow==3&&ptHigh==4.5) ctauCut=0.1135;
  //  else if(ptLow==4.5&&ptHigh==6.5) ctauCut=0.0935;
  //  else if(ptLow==6.5&&ptHigh==9) ctauCut=0.0755;
  //  else if(ptLow==9&&ptHigh==12)  ctauCut=0.0635;
  //  else if(ptLow==12&&ptHigh==30) ctauCut=0.0495;}
  cout<<"pt["<<ptLow<<" - "<<ptHigh<<" GeV/c], "<<"ctau cut: "<<ctauCut<<endl;

  const int nMassBin = 12;

  float massBinDiff[nMassBin+1]={2.6, 2.7, 2.8, 2.9, 3.0, 3.06, 3.09, 3.12, 3.15, 3.2, 3.3, 3.4, 3.5};
  //float massBinDiff[nMassBin+1]={2.6, 2.75, 2.9, 3.0, 3.06, 3.09, 3.12, 3.15, 3.2,3.5};
  float massBin_[nMassBin+1];

  kineLabel = kineLabel + Form("_m%.1f-%.1f",massLow,massHigh) + "_" + dimusignString;
  cout<<kineLabel<<endl;

  for(int i=0; i<=nMassBin; i++){
    massBin_[i] = massBinDiff[i];
    //massBin[i] = 2.6 + massBinDiff[i]*0.02;
  }

  float* massBin;
  int nMassBin_;
  massBin = massBin_; nMassBin_ = nMassBin;

  int bfevt =-1;
  int afevt =-1;

  double mass_low_SB1 = 2.6; 
  double mass_high_SB1 = 2.9; 
  double mass_low_SB = 2.9;//3.0 
  double mass_high_SB = 3.3; 
  double mass_low_SB2 = 3.3; //3.3
  double mass_high_SB2 = 3.5; 

  int nQBin = 200;
  const int nMass_div = 3;

  TString fSB[nMass_div] = {"SB1 (2.6<m<2.9)","SB2 (3.3<m<3.5)","S (2.9<m<3.3)"};

  //Define drawing histogram
  TH1D* h_v2_1[nMass_div];
  TH1D* h_v2_2[nMass_div];
  TH1D* h_v2_3[nMass_div];
  TH1D* h_v2_4[nMass_div];
  TH1D* h_ljpsi[nMass_div];

  double Q_avg_low = -6500;
  double Q_avg_high = 27000;
  double Q_avg_low_dimu = -150;
  double Q_avg_high_dimu = 150;
  for(int imass=0; imass<nMass_div;imass++){
    h_v2_1[imass] = new TH1D(Form("h_v2_1_%d",imass),";#LTQ_{2}Q_{2A}^{*}#GT;counts",nQBin,Q_avg_low_dimu,Q_avg_high_dimu);
    h_v2_2[imass] = new TH1D(Form("h_v2_2_%d",imass),";#LTQ_{2A}Q_{2B}^{*}#GT;counts",nQBin,Q_avg_low,Q_avg_high);
    h_v2_3[imass] = new TH1D(Form("h_v2_3_%d",imass),";#LTQ_{2A}Q_{2C}^{*}#GT;counts",nQBin,Q_avg_low,Q_avg_high);
    h_v2_4[imass] = new TH1D(Form("h_v2_4_%d",imass),";#LTQ_{2B}Q_{2C}^{*}#GT;counts",nQBin,Q_avg_low,Q_avg_high);
    h_ljpsi[imass] = new TH1D(Form("h_ljBkg_%d",imass),";l_{J/#psi(mm)};Counts/(0.03 mm)",1050,-1.5,2);
  }

  TH1D* h_v2_num_q1 = new TH1D("h_v2_num_q1",";m_{#mu^{+}#mu^{-}};#LTQ_{2}Q_{2A}^{*}#GT" ,nMassBin,massBin);
  TH1D* h_v2_den_q2 = new TH1D("h_v2_num_q2",";m_{#mu^{+}#mu^{-}};#LTQ_{2A}Q_{2B}^{*}#GT",nMassBin,massBin);
  TH1D* h_v2_den_q3 = new TH1D("h_v2_num_q3",";m_{#mu^{+}#mu^{-}};#LTQ_{2A}Q_{2C}^{*}#GT",nMassBin,massBin);
  TH1D* h_v2_den_q4 = new TH1D("h_v2_num_q4",";m_{#mu^{+}#mu^{-}};#LTQ_{2B}Q_{2C}^{*}#GT",nMassBin,massBin);

  TH1D* h_mass = new TH1D("h_mass",";m_{#mu^{+}#mu^{-}};Counts/(0.02 GeV)",120,massLow,massHigh);
  TH1D* h_decayL = new TH1D("h_decayL",";l_{J/#psi(mm)};Counts/(0.03 mm)",4000,-1,3);
  //TH1D* h_ljpsi_bkg = new TH1D("h_ljpsi_bkg",";l_{J/#psi(mm)};Counts/(0.03 mm)",1050,-1.5,2);

  RooRealVar* massVar = new RooRealVar("mass","mass variable",0,200,"GeV/c^{2}");
  RooRealVar* weightVar = new RooRealVar("weight", "weight var", 0, 10000,"");

  RooArgSet* argSet = new RooArgSet(*massVar, *weightVar);
  RooDataSet* dataSet = new RooDataSet("dataset","a dataset",*argSet); 

  const static int countMax = 1000000;
  vector<vector<double>> v2_1(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v2_2(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v2_3(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v2_4(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v2_1_raw(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v2_2_raw(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v2_3_raw(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v2_4_raw(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> weight_dimu(nMassBin,vector<double> (countMax,0));

  vector<double> v2_1_avg(nMassBin,0);
  vector<double> v2_2_avg(nMassBin,0);
  vector<double> v2_3_avg(nMassBin,0);
  vector<double> v2_4_avg(nMassBin,0);

  int dbcount=0;
  vector<unsigned int> count(nMassBin,0);
  vector<unsigned int> count_ss(nMassBin,0);

  vector<double> weight_s(nMassBin,0);

  double weight_acc = 1;
  double weight_eff = 1;

  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;

  //Begin Loop
  for(int i=0; i<nEvt; i++){
    tree->GetEntry(i);
    if(!(cBin>=cLow||cBin<=cHigh)) continue;
    if(fabs(vz)>=15) continue;
    // Fill Dimuon Loop
    for(int j=0; j<nDimu; j++){
      if(pt[j]>ptLow&&pt[j]<ptHigh&&recoQQsign[j]==0&&mass[j]>=massLow&&mass[j]<=massHigh&&abs(y[j])>yLow&&abs(y[j])<yHigh&&
          //abs(eta1[j])>1.6&&abs(eta1[j])<2.4&&abs(eta2[j])>1.6&&abs(eta2[j])<2.4&&
			  ( ((abs(eta1[j]) <= 1.2) && (pt1[j] >=3.5)) || ((abs(eta2[j]) <= 1.2) && (pt2[j] >=3.5)) ||
				((abs(eta1[j]) > 1.2) && (abs(eta1[j]) <= 2.1) && (pt1[j] >= 5.47-1.89*(abs(eta1[j])))) || ((abs(eta2[j]) > 1.2)  && (abs(eta2[j]) <= 2.1) && (pt2[j] >= 5.47-1.89*(abs(eta2[j])))) ||
				((abs(eta1[j]) > 2.1) && (abs(eta1[j]) <= 2.4) && (pt1[j] >= 1.5)) || ((abs(eta2[j]) > 2.1)  && (abs(eta2[j]) <= 2.4) && (pt2[j] >= 1.5)) )
		){
		  //        if ((PR==0||PR==1) && ctau3D[j]>ctauCut) continue;
		  weight_acc=1;
		  weight_eff=1;
		  if(fAccW){
			  if(cBin<20) {
				  if ( abs((double)y[j])<0.6) { weight_acc = getAccWeight(hAccPt[0], pt[j]); }
				  else if ( abs((double)y[j])>=0.6 && abs((double)y[j])<1.2 ) {weight_acc = getAccWeight(hAccPt[1], pt[j]); }
				  else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) {weight_acc = getAccWeight(hAccPt[2], pt[j]); }
				  else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<2.4 ) {weight_acc = getAccWeight(hAccPt[3], pt[j]); }
			  }
			  if(cBin>=20 && cBin<100){
				  if ( abs((double)y[j])<0.6) { weight_acc = getAccWeight(hAccPt[4], pt[j]); }
				  else if ( abs((double)y[j])>=0.6 && abs((double)y[j])<1.2 ) {weight_acc = getAccWeight(hAccPt[5], pt[j]); }
				  else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) {weight_acc = getAccWeight(hAccPt[6], pt[j]); }
				  else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<2.4 ) {weight_acc = getAccWeight(hAccPt[7], pt[j]); }
			  }
			  if(cBin>=100 && cBin<180){
				  if ( abs((double)y[j])<0.6) { weight_acc = getAccWeight(hAccPt[8], pt[j]); }
				  else if ( abs((double)y[j])>=0.6 && abs((double)y[j])<1.2 ) {weight_acc = getAccWeight(hAccPt[9], pt[j]); }
				  else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) {weight_acc = getAccWeight(hAccPt[10], pt[j]); }
				  else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<2.4 ) {weight_acc = getAccWeight(hAccPt[11], pt[j]); }
			  }
		  }
		  if(fEffW){
			  if(cBin<20) {
				  if ( abs((double)y[j])<0.6) { weight_eff = getEffWeight(hEffPt[0], pt[j]); }
				  else if ( abs((double)y[j])>=0.6 && abs((double)y[j])<1.2 ) {weight_eff = getEffWeight(hEffPt[1], pt[j]); }
				  else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) {weight_eff = getEffWeight(hEffPt[2], pt[j]); }
				  else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<2.4 ) {weight_eff = getEffWeight(hEffPt[3], pt[j]); }
			  }
			  if(cBin>=20 && cBin<100){
				  if ( abs((double)y[j])<0.6) { weight_eff = getEffWeight(hEffPt[4], pt[j]); }
				  else if ( abs((double)y[j])>=0.6 && abs((double)y[j])<1.2 ) {weight_eff = getEffWeight(hEffPt[5], pt[j]); }
				  else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) {weight_eff = getEffWeight(hEffPt[6], pt[j]); }
				  else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<2.4 ) {weight_eff = getEffWeight(hEffPt[7], pt[j]); }
			  }
			  if(cBin>=100 && cBin<180){
				  if ( abs((double)y[j])<0.6) { weight_eff = getEffWeight(hEffPt[8], pt[j]); }
				  else if ( abs((double)y[j])>=0.6 && abs((double)y[j])<1.2 ) {weight_eff = getEffWeight(hEffPt[9], pt[j]); }
				  else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) {weight_eff = getEffWeight(hEffPt[10], pt[j]); }
				  else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<2.4 ) {weight_eff = getEffWeight(hEffPt[11], pt[j]); }
			  }
		  }
		double weight_ = weight * weight_eff * weight_acc;
        for(int imbin=0; imbin<nMassBin; imbin++){
          if(mass[j]>=massBin[imbin] && mass[j]<massBin[imbin+1]){
            v2_1[imbin][count[imbin]] = (qxa[j]*qxdimu[j] + qya[j]*qydimu[j])*weight_;
            v2_2[imbin][count[imbin]] = (qxa[j]*qxb[j] + qya[j]*qyb[j])*weight_;
            v2_3[imbin][count[imbin]] = (qxa[j]*qxc[j] + qya[j]*qyc[j])*weight_;
            v2_4[imbin][count[imbin]] = (qxb[j]*qxc[j] + qyb[j]*qyc[j])*weight_;
            v2_1_raw[imbin][count[imbin]] = (qxa[j]*qxdimu[j] + qya[j]*qydimu[j]);
            v2_2_raw[imbin][count[imbin]] = (qxa[j]*qxb[j] + qya[j]*qyb[j]);
            v2_3_raw[imbin][count[imbin]] = (qxa[j]*qxc[j] + qya[j]*qyc[j]);
            v2_4_raw[imbin][count[imbin]] = (qxb[j]*qxc[j] + qyb[j]*qyc[j]);

            v2_1_avg[imbin] += v2_1[imbin][count[imbin]];
            v2_2_avg[imbin] += v2_2[imbin][count[imbin]];
            v2_3_avg[imbin] += v2_3[imbin][count[imbin]];
            v2_4_avg[imbin] += v2_4[imbin][count[imbin]];

			weight_dimu[imbin][count[imbin]] = weight_;

            count[imbin]++;
			weight_s[imbin] += weight_;
          }
        }
        if(mass[j]>=mass_low_SB1 && mass[j]<mass_high_SB1){
          h_v2_1[0]->Fill(qxa[j]*qxdimu[j] + qya[j]*qydimu[j], weight_);
          h_v2_2[0]->Fill(qxa[j]*qxb[j] + qya[j]*qyb[j], weight_);
          h_v2_3[0]->Fill(qxa[j]*qxc[j] + qya[j]*qyc[j], weight_);
          h_v2_4[0]->Fill(qxb[j]*qxc[j] + qyb[j]*qyc[j], weight_);
          h_ljpsi[0]->Fill(ctau3D[j]);
        }
        else if(mass[j]>=mass_low_SB2 && mass[j]<mass_high_SB2){
          h_v2_1[1]->Fill(qxa[j]*qxdimu[j] + qya[j]*qydimu[j], weight_);
          h_v2_2[1]->Fill(qxa[j]*qxb[j] + qya[j]*qyb[j], weight_);
          h_v2_3[1]->Fill(qxa[j]*qxc[j] + qya[j]*qyc[j], weight_);
          h_v2_4[1]->Fill(qxb[j]*qxc[j] + qyb[j]*qyc[j], weight_);
          h_ljpsi[1]->Fill(ctau3D[j]);
        }
        else if(mass[j]>=mass_low_SB && mass[j]<mass_high_SB){
          h_v2_1[2]->Fill(qxa[j]*qxdimu[j] + qya[j]*qydimu[j], weight_);
          h_v2_2[2]->Fill(qxa[j]*qxb[j] + qya[j]*qyb[j], weight_);
          h_v2_3[2]->Fill(qxa[j]*qxc[j] + qya[j]*qyc[j], weight_);
          h_v2_4[2]->Fill(qxb[j]*qxc[j] + qyb[j]*qyc[j], weight_);
          h_ljpsi[2]->Fill(ctau3D[j]);
        }
        h_decayL->Fill(ctau3D[j]);
        massVar->setVal((double)mass[j]);
		weightVar->setVal((double)weight_);
        dataSet->add(*argSet);
        h_mass->Fill(mass[j],weight_);
      }
    }
  }
  cout<<"How many Jpsi??: "<<float(h_mass->GetEntries())<<endl;

  TLine *lcut;
  //for(int bins=0; bins<h_decayL->GetNbinsX(); bins++){
  //  h_decayL->Integral(0,bins);
  //  if(h_decayL->Integral(0,bins)>=0.79*h_decayL->GetEntries()
  //      &&h_decayL->Integral(0,bins)<=0.81*h_decayL->GetEntries()){ 
  //    cout<<"THE BIN IS: "<<bins<<endl;
  //    lcut = new TLine(h_decayL->GetBinCenter(bins), 0, h_decayL->GetBinCenter(bins), h_decayL->GetMaximum());
  //  }
  //}
  //for(int bins=0; bins<h_decayL->GetNbinsX(); bins++){
  //  h_decayL->Integral(0,bins);
  //  if(h_decayL->Integral(0,bins)<=0.89*h_decayL->GetEntries() || h_decayL->Integral(0,bins)>0.90*h_decayL->GetEntries()) continue;
  //    cout<<"ANSWER: "<<bins<<"th bin"<<", ljpsi: "<<h_decayL->GetBinCenter(bins)<<endl;
  //    lcut = new TLine(h_decayL->GetBinCenter(bins), 0, h_decayL->GetBinCenter(bins), 1);
  //}

  int nMassFrameBin = 120;
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataSet);
  ws->data("dataset")->Print();
  RooDataSet *reducedDS = new RooDataSet("reducedDS","A sample",*dataSet->get(),Import(*dataSet),WeightVar(*ws->var("weight")));
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();

  RooPlot* myPlot = ws->var("mass")->frame(nMassFrameBin);
  ws->data("reducedDS")->plotOn(myPlot,Name("massDataHist"));

  RooHist* hist = (RooHist*) myPlot->findObject("massDataHist");
  int nHistP = hist->GetN();
  cout << "nHistP : " << nHistP << endl;
  if(nHistP != nMassFrameBin) cout << "ERROR::: INCONSISTENT NUMBER OF BINS" << endl;

  TGraphAsymmErrors *g_mass = new TGraphAsymmErrors();
  g_mass->SetName("g_mass");
  Double_t xp, yp, ypl, yph;
  for(int ip = 0; ip < nHistP; ip++){
    xp=0; yp=0; ypl=0; yph=0;
    hist->GetPoint(ip,xp,yp);
    ypl = hist->GetErrorYlow(ip);
    yph = hist->GetErrorYhigh(ip);
    g_mass->SetPoint(ip,xp,yp);
    g_mass->SetPointError(ip,0,0,ypl,yph);
  }

  TF1* fitmc1;

  vector<double> v2_1_err(nMassBin,0);
  vector<double> v2_2_err(nMassBin,0);
  vector<double> v2_3_err(nMassBin,0);
  vector<double> v2_4_err(nMassBin,0);

  for(int ibin=0; ibin<nMassBin; ibin++){
    /*v2_1_avg[ibin] = v2_1_avg[ibin]/count[ibin];
    v2_2_avg[ibin] = v2_2_avg[ibin]/count[ibin];
    v2_3_avg[ibin] = v2_3_avg[ibin]/count[ibin];
    v2_4_avg[ibin] = v2_4_avg[ibin]/count[ibin];
	*/
    v2_1_avg[ibin] = v2_1_avg[ibin]/weight_s[ibin];
    v2_2_avg[ibin] = v2_2_avg[ibin]/weight_s[ibin];
    v2_3_avg[ibin] = v2_3_avg[ibin]/weight_s[ibin];
    v2_4_avg[ibin] = v2_4_avg[ibin]/weight_s[ibin];

    for(int icount=0; icount<count[ibin]; icount++){
      v2_1_err[ibin] += (v2_1_raw[ibin][icount]-v2_1_avg[ibin])*(v2_1_raw[ibin][icount]-v2_1_avg[ibin]) * weight_dimu[ibin][icount] * weight_dimu[ibin][icount];
      v2_2_err[ibin] += (v2_2_raw[ibin][icount]-v2_2_avg[ibin])*(v2_2_raw[ibin][icount]-v2_2_avg[ibin]) * weight_dimu[ibin][icount] * weight_dimu[ibin][icount];
      v2_3_err[ibin] += (v2_3_raw[ibin][icount]-v2_3_avg[ibin])*(v2_3_raw[ibin][icount]-v2_3_avg[ibin]) * weight_dimu[ibin][icount] * weight_dimu[ibin][icount];
      v2_4_err[ibin] += (v2_4_raw[ibin][icount]-v2_4_avg[ibin])*(v2_4_raw[ibin][icount]-v2_4_avg[ibin]) * weight_dimu[ibin][icount] * weight_dimu[ibin][icount];
    }

    /*v2_1_err[ibin] = TMath::Sqrt(v2_1_err[ibin])/count[ibin];
    v2_2_err[ibin] = TMath::Sqrt(v2_2_err[ibin])/count[ibin];
    v2_3_err[ibin] = TMath::Sqrt(v2_3_err[ibin])/count[ibin];
    v2_4_err[ibin] = TMath::Sqrt(v2_4_err[ibin])/count[ibin];
	*/
    v2_1_err[ibin] = TMath::Sqrt(v2_1_err[ibin])/weight_s[ibin];
    v2_2_err[ibin] = TMath::Sqrt(v2_2_err[ibin])/weight_s[ibin];
    v2_3_err[ibin] = TMath::Sqrt(v2_3_err[ibin])/weight_s[ibin];
    v2_4_err[ibin] = TMath::Sqrt(v2_4_err[ibin])/weight_s[ibin];

    h_v2_num_q1->SetBinContent(ibin+1,v2_1_avg[ibin]);
    h_v2_num_q1->SetBinError(ibin+1,v2_1_err[ibin]);
    h_v2_den_q2->SetBinContent(ibin+1,v2_2_avg[ibin]);
    h_v2_den_q2->SetBinError(ibin+1,v2_2_err[ibin]);
    h_v2_den_q3->SetBinContent(ibin+1,v2_3_avg[ibin]);
    h_v2_den_q3->SetBinError(ibin+1,v2_3_err[ibin]);
    h_v2_den_q4->SetBinContent(ibin+1,v2_4_avg[ibin]);
    h_v2_den_q4->SetBinError(ibin+1,v2_4_err[ibin]);

    cout << ibin << "th Bin : " << count[ibin] << ",  weight sum  : " << weight_s[ibin] << endl;
    cout << "v2_1_avg " << ibin << " : " << v2_1_avg[ibin] << endl;
    cout << "v2_2_avg " << ibin << " : " << v2_2_avg[ibin] << endl;
    cout << "v2_3_avg " << ibin << " : " << v2_3_avg[ibin] << endl;
    cout << "v2_4_avg " << ibin << " : " << v2_4_avg[ibin] << endl;

    cout << "h_v2_num_q1 " << ibin << ", val : " << h_v2_num_q1->GetBinContent(ibin+1) << " err : " << h_v2_num_q1->GetBinError(ibin+1) << endl;
    cout << "h_v2_den_q2 " << ibin << ", val : " << h_v2_den_q2->GetBinContent(ibin+1) << " err : " << h_v2_den_q2->GetBinError(ibin+1) << endl;
    cout << "h_v2_den_q3 " << ibin << ", val : " << h_v2_den_q3->GetBinContent(ibin+1) << " err : " << h_v2_den_q3->GetBinError(ibin+1) << endl;
    cout << "h_v2_den_q4 " << ibin << ", val : " << h_v2_den_q4->GetBinContent(ibin+1) << " err : " << h_v2_den_q4->GetBinError(ibin+1) << endl;
  }


  TH1D* h_v2_den_ = (TH1D*)h_v2_den_q2->Clone("h_v2_den_");
  h_v2_den_->Multiply(h_v2_den_q3);
  h_v2_den_->Divide(h_v2_den_q4);

  TH1D* h_v2_den = (TH1D*) h_v2_den_->Clone("h_v2_den"); //h_v2_den->Reset();
  GetHistSqrt(h_v2_den_,h_v2_den);

  TH1D* h_v2_final = (TH1D*) h_v2_num_q1 -> Clone("h_v2_SplusB");
  h_v2_final->Divide(h_v2_den);

  SetHistStyle(h_v2_final,0,0);
  SetGraphStyle2(g_mass,0,0);
  //g_mass->GetYaxis()->SetLimits(0,10000);
  //g_mass->GetYaxis()->SetRangeUser(0,);
  g_mass->GetXaxis()->SetLimits(massLow,massHigh);
  g_mass->GetXaxis()->SetRangeUser(massLow,massHigh);
  g_mass->GetYaxis()->SetTitle("Events/(0.02 GeV)");
  g_mass->SetMarkerSize(1);
  g_mass->GetYaxis()->SetLabelSize(0.055);
  g_mass->GetXaxis()->SetLabelSize(0.04);
  g_mass->GetYaxis()->SetTitleSize(0.055);
  g_mass->GetYaxis()->SetTitleOffset(1.7);
  g_mass->GetXaxis()->SetTitleSize(0.065);
  h_v2_final->GetYaxis()->SetRangeUser(-0.05,0.26);
  h_v2_final->GetYaxis()->SetTitle("v_{2}^{S+B}");
  h_v2_final->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
  h_v2_final->GetYaxis()->SetLabelSize(0.055);
  h_v2_final->GetXaxis()->SetLabelSize(0.055);
  h_v2_final->GetXaxis()->SetTitleSize(0.07);
  h_v2_final->GetYaxis()->SetTitleSize(0.07);
  h_v2_final->GetYaxis()->SetTitleOffset(1.2);
  h_v2_final->GetXaxis()->SetTitleOffset(1.0);
  h_v2_final->GetXaxis()->SetLabelOffset(0.011);
  TLine *l=new TLine(2.6, 0, 3.5, 0);
  l->SetLineColor(kBlack);

  float pos_x = 0.43;
  float pos_x_mass = 0.73;
  float pos_y = 0.75;
  float pos_y_diff = 0.071;
  int text_color = 1;
  float text_size = 16;
  TString perc = "%";

  //TGaxis::SetMaxDigits(2);
  TCanvas* c_mass_v2 = new TCanvas("c_mass_v2","",590,750);
  TPad* pad1 = new TPad("pad1","pad1",0,0.5,1.0,1.0);
  TPad* pad2 = new TPad("pad2","pad2",0,0.0,1.0,0.5);
  c_mass_v2->cd();
  pad1->SetTicks(1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLeftMargin(0.19);
  pad1->SetTopMargin(0.08);
  pad1->cd();
  g_mass->GetXaxis()->SetTitleSize(0);
  g_mass->GetXaxis()->SetLabelSize(0);
  g_mass->Draw("AP");
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  else if(ptLow ==6.5) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  else if(ptHigh==6.5) drawText(Form("%.f < p_{T}^{#mu#mu} < %.1f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  else if(ptLow!=0) drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
  drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);
  //    drawText(Form("L_{J/psi} < %.4f", ctauCut),pos_x_mass,pos_y-pos_y_diff*5,text_color,text_size);
  if(PR==0){drawText(Form("L_{J/psi} %.4f", ctauCut),pos_x_mass,pos_y-pos_y_diff*5,text_color,text_size);}
  else if(PR==1){drawText(Form("L_{J/psi} %.4f", ctauCut),pos_x_mass,pos_y-pos_y_diff*5,text_color,text_size);}
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.17);
  pad2->SetLeftMargin(0.19);
  pad2->SetTicks();
  pad2->cd();
  ////////////pad2//////////
  jumSun(massLow,0,massHigh,0,1,1);
  h_v2_final->Draw("P");
  l->Draw("same");

  //CMS_lumi_v2mass(pad1,iPeriod,iPos,1);  
  CMS_lumi_v2mass(pad1,iPeriod,iPos);  
  pad1->Update();
  pad2->Update();
  c_mass_v2->cd();
  pad1->Draw();
  pad2->Draw();
  //gSystem->mkdir(Form("figs/q_vector",sample.Data()),1);
  c_mass_v2->SaveAs(Form("figs/v2mass_hist/v2Mass_%s_%s.pdf", bCont.Data(), kineLabel.Data()));

  TCanvas* c_decayL = new TCanvas("c_decayL","",600,600);
  c_decayL->SetLogy();
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  else if(ptLow!=0) drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
  drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);
  h_decayL->Draw("P");
  //lcut->Draw("same");
  c_decayL->Update();
  c_decayL->SaveAs(Form("figs/decayL/decayL_%s_%s.pdf", bCont.Data(), kineLabel.Data()));

  TCanvas* c_mass = new TCanvas("c_mass","",600,600);
  c_mass->cd();
  h_mass->Draw("P");
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  else if(ptLow!=0) drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
  drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);
  //   drawText(Form("L_{J/psi} < %.4f", ctauCut),pos_x_mass,pos_y-pos_y_diff*5,text_color,text_size);
  drawText(Form("L_{J/psi} > %.4f", ctauCut),pos_x_mass,pos_y-pos_y_diff*5,text_color,text_size);
  //CMS_lumi_v2mass(c_mass,iPeriod,iPos,0);
  CMS_lumi_v2mass(c_mass,iPeriod,iPos);
  c_mass->Update();
  c_mass->SaveAs(Form("figs/mass_dist/MassDist_%s_%s.pdf", bCont.Data(), kineLabel.Data()));

  TLegend *leg_v2_1 = new TLegend(0.38,0.64,0.77,0.9);
  TLegend *leg_v2_2 = new TLegend(0.38,0.64,0.77,0.9);
  TLegend *leg_v2_3 = new TLegend(0.38,0.64,0.77,0.9);
  TLegend *leg_v2_4 = new TLegend(0.38,0.64,0.77,0.9);
  SetLegendStyle(leg_v2_1);
  SetLegendStyle(leg_v2_2);
  SetLegendStyle(leg_v2_3);
  SetLegendStyle(leg_v2_4);

  double mean_v2_1[nMass_div];
  double mean_v2_2[nMass_div];
  double mean_v2_3[nMass_div];
  double mean_v2_4[nMass_div];

  double mean_err_v2_1[nMass_div];
  double mean_err_v2_2[nMass_div];
  double mean_err_v2_3[nMass_div];
  double mean_err_v2_4[nMass_div];

  for(int i=0;i<nMass_div;i++){
    SetHistStyle(h_v2_1[i],i,i);
    SetHistStyle(h_v2_2[i],i,i);
    SetHistStyle(h_v2_3[i],i,i);
    SetHistStyle(h_v2_4[i],i,i);

    scaleInt(h_v2_1[i]);
    scaleInt(h_v2_2[i]);
    scaleInt(h_v2_3[i]);
    scaleInt(h_v2_4[i]);

    mean_v2_1[i] = h_v2_1[i]->GetMean();
    mean_v2_2[i] = h_v2_2[i]->GetMean();
    mean_v2_3[i] = h_v2_3[i]->GetMean();
    mean_v2_4[i] = h_v2_4[i]->GetMean();

    mean_err_v2_1[i] = h_v2_1[i]->GetMeanError();
    mean_err_v2_2[i] = h_v2_2[i]->GetMeanError();
    mean_err_v2_3[i] = h_v2_3[i]->GetMeanError();
    mean_err_v2_4[i] = h_v2_4[i]->GetMeanError();

    leg_v2_1->AddEntry(h_v2_1[i],(fSB[i]+Form(" mean:%.2f",mean_v2_1[i])+Form(" err:%.2f",mean_err_v2_1[i])).Data(),"l");
    leg_v2_2->AddEntry(h_v2_2[i],(fSB[i]+Form(" mean:%.2f",mean_v2_2[i])+Form(" err:%.2f",mean_err_v2_2[i])).Data(),"l");
    leg_v2_3->AddEntry(h_v2_3[i],(fSB[i]+Form(" mean:%.2f",mean_v2_3[i])+Form(" err:%.2f",mean_err_v2_3[i])).Data(),"l");
    leg_v2_4->AddEntry(h_v2_4[i],(fSB[i]+Form(" mean:%.2f",mean_v2_4[i])+Form(" err:%.2f",mean_err_v2_4[i])).Data(),"l");
  }
  //cout<<"what?? : "<<h_v2_1[0]->GetEntries()<<endl;

  //TCanvas *c_ljBkg = new TCanvas("c_ljBkg","",600,600);
  //c_ljBkg->cd();
  //h_ljpsi[0]->Draw("hist");
  //h_ljpsi[0]->SetLineColor(kBlue+2);
  //h_ljpsi[1]->Draw("hist same");
  //h_ljpsi[1]->SetLineColor(kRed+2);
  //h_ljpsi[2]->Draw("hist same");
  //h_ljpsi[2]->SetLineColor(kBlack);

  //TH1D* h_ljpsi_bkg = (TH1D*) h_ljpsi[0]->Clone("h_ljpsi_bkg");
  //h_ljpsi_bkg->Scale(h_ljpsi[0]->GetEntries()/(h_ljpsi[0]->GetEntries()+h_ljpsi[1]->GetEntries()));
  //h_ljpsi_bkg->Draw("EP same");
  //h_ljpsi_bkg->SetMarkerColor(kRed+2);
  //h_ljpsi_bkg->SetMarkerStyle(4);

  TCanvas *c_qq_1 = new TCanvas("c_qqa","",600,600);
  c_qq_1->cd();
  h_v2_1[2]->Draw("hist");
  h_v2_1[1]->Draw("hist same");
  h_v2_1[0]->Draw("hist same");
  leg_v2_1->Draw("same");
  leg_v2_1->SetFillStyle(0);
  c_qq_1->SaveAs(Form("figs/q_vector/c_qqa_%s_%s_Eff%d_Acc%d_PtW%d_TnP%d.pdf",bCont.Data(), kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));

  TCanvas *c_qq_2 = new TCanvas("c_qaqb","",600,600);
  c_qq_2->cd();
  h_v2_2[2]->Draw("hist");
  h_v2_2[1]->Draw("hist same");
  h_v2_2[0]->Draw("hist same");
  leg_v2_2->Draw("same");
  c_qq_2->SaveAs(Form("figs/q_vector/c_qaqb_%s_%s_Eff%d_Acc%d_PtW%d_TnP%d.pdf",bCont.Data(), kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));

  TCanvas *c_qq_3 = new TCanvas("c_qaqc","",600,600);
  c_qq_3->cd();
  h_v2_3[2]->Draw("hist");
  h_v2_3[1]->Draw("hist same");
  h_v2_3[0]->Draw("hist same");
  leg_v2_3->Draw("same");
  c_qq_3->SaveAs(Form("figs/q_vector/c_qaqc_%s_%s_Eff%d_Acc%d_PtW%d_TnP%d.pdf",bCont.Data(), kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));

  TCanvas *c_qq_4 = new TCanvas("c_qbqc","",600,600);
  c_qq_4->cd();
  h_v2_4[2]->Draw("hist");
  h_v2_4[1]->Draw("hist same");
  h_v2_4[0]->Draw("hist same");
  leg_v2_4->Draw("same");
  c_qq_4->SaveAs(Form("figs/q_vector/c_qbqc_%s_%s_Eff%d_Acc%d_PtW%d_TnP%d.pdf",bCont.Data(), kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));

  TFile *wf = new TFile(Form("roots/v2mass_hist/Jpsi_%s_%s_Eff%d_Acc%d_PtW%d_TnP%d.root", bCont.Data(), kineLabel.Data(),fEffW,fAccW,isPtW,isTnP),"recreate");
  wf->cd();
  h_v2_final->Write();
  h_decayL->Write();
  g_mass->Write();
  h_mass->Write();
}

void GetHistSqrt(TH1D* h1, TH1D* h2){
  if(h1->GetNbinsX() != h2->GetNbinsX()){ cout << "Inconsistent # of bins b/w histograms !! " << endl;}
  double content;
  double err;
  for(int i=1; i<=h1->GetNbinsX(); i++){
    content=0;err=0;
    content = h1->GetBinContent(i);
    err = h1->GetBinError(i);
    err = 0.5*err*TMath::Power(content,-0.5);
    h2->SetBinContent(i,TMath::Sqrt(content));
    h2->SetBinError(i,err);
  }
} 

void GetHistBkg(TH1D* h1, TH1D* h2){
  double Nh1;
  double Nh2;
  double Nh12;
  Nh1 = h1->GetEntries();
  Nh2 = h2->GetEntries();
  Nh12 = Nh1+Nh2;
  return Nh1/(Nh1+Nh2);

}

double getAccWeight(TH1D* h, double pt){
  double binN = h->FindBin(pt);
  TF1 *acc1 = (TF1*)h->FindObject("f1");
  double acc = acc1->Eval(pt);
  double weight_ = 1./acc;
  return weight_;
}

double getEffWeight(TH1D *h, double pt){
  double binN = h->FindBin(pt);
  TF1 *eff1 = (TF1*)h->FindObject("f1");
  double eff = eff1->Eval(pt);
  double weight_ = 1./eff;
  return weight_;
}
