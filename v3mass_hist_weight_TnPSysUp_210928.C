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

void v3mass_hist_weight_TnPSysUp_210928(
    double ptLow =  6.5, double ptHigh =  7.5,
    double yLow = 0, double yHigh = 2.4,
    int cLow = 20, int cHigh = 120,
    int ctauCut = 2, //-1:Left, 0:Central, 1:Right 2:w/o cut
    float massLow = 2.6, float massHigh =3.5, 
    bool dimusign=true, bool isMC = false, 
    bool fAccW = true, bool fEffW = true,
    bool isTnP = true, bool isPtW = true,
    int dtype = 1, 
    int weight_PR = 0, // PR : 0, NP : 1
    int PR=2
    ) //PR 0: PR, 1: NP, 2: Inc.
{
  //Basic Setting
  gStyle->SetOptStat(0);
  //TString DATE="210507";
  //TString DATE="Corr";
  //TString DATE="210503";
  //TString DATE="10_60";
  //TString DATE="20_40";
  //TString DATE="0_180";
  //TString DATE="210520";
  TString DATE;
  if(ptLow==6.5&&ptHigh==50) DATE=Form("%i_%i",0,180);
  else DATE=Form("%i_%i",cLow/2,cHigh/2);
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir(Form("figs/TnPSysUp_210928/v3mass_hist/%s",DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/TnPSysUp_210928/q_vector/%s",DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/TnPSysUp_210928/mass_dist/%s",DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/TnPSysUp_210928/decayL/%s",DATE.Data()), kTRUE);
  gSystem->mkdir(Form("roots/TnPSysUp_210928/v3mass_hist/%s",DATE.Data()), kTRUE);
  gStyle->SetOptStat(000000000);
  gROOT->ForceStyle();
  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  TH1::SetDefaultSumw2();
  //TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh, 1) ;
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0, cLow, cHigh) ;
  TString dimusignString;
  if(dimusign) dimusignString = "OS";
  else if(!dimusign) dimusignString = "SS";
  TString sample;
  TString wName;
  TString cutName;

  if(weight_PR==0) wName="prompt";
  else if(weight_PR==1) wName="nprompt";

  cout<<"Data type: "<<dtype<<endl;

  TChain *tree = new TChain("mmepevt");
  if(!isMC){
    TString f1 = "skimmedFiles/OniaFlowSkim_JpsiTrig_DBAllPD_isMC0_HFNom_v3_210531.root";
    //TString f1 = "skimmedFiles/OniaFlowSkim_JpsiTrig_AllDBPD_isMC0_HFNom_210726.root";
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

  TFile *fEff;
  if(ptLow==6.5&&ptHigh==50) fEff= new TFile(Form("./Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_%s_pbpb_Jpsi_PtW%d_tnp%d_210913_TnPSys.root",wName.Data(),isPtW,isTnP),"read");
  else fEff= new TFile(Form("./Eff_Acc/roots/mc_eff_vs_pt_cent_%i_to_%i_rap_%s_pbpb_Jpsi_PtW%d_tnp%d_210913_TnPSys.root",cLow,cHigh,wName.Data(),isPtW,isTnP),"read");

  TH1D* hEffPt[4];
  if(ptLow==6.5&&ptHigh==50){
    hEffPt[0] = (TH1D*) fEff -> Get(Form("hEff_sysUp_mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_180_absy0_1p2",isTnP,isPtW));
    hEffPt[1] = (TH1D*) fEff -> Get(Form("hEff_sysUp_mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_180_absy1p2_1p6",isTnP,isPtW));
    hEffPt[2] = (TH1D*) fEff -> Get(Form("hEff_sysUp_mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_180_absy1p6_1p8",isTnP,isPtW));
    hEffPt[3] = (TH1D*) fEff -> Get(Form("hEff_sysUp_mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_180_absy1p8_2p4",isTnP,isPtW));}
  else{
    hEffPt[0] = (TH1D*) fEff -> Get(Form("hEff_sysUp_mc_eff_vs_pt_TnP%d_PtW%d_cent_%i_to_%i_absy0_1p2",isTnP,isPtW,cLow,cHigh));
    hEffPt[1] = (TH1D*) fEff -> Get(Form("hEff_sysUp_mc_eff_vs_pt_TnP%d_PtW%d_cent_%i_to_%i_absy1p2_1p6",isTnP,isPtW,cLow,cHigh));
    hEffPt[2] = (TH1D*) fEff -> Get(Form("hEff_sysUp_mc_eff_vs_pt_TnP%d_PtW%d_cent_%i_to_%i_absy1p6_1p8",isTnP,isPtW,cLow,cHigh));
    hEffPt[3] = (TH1D*) fEff -> Get(Form("hEff_sysUp_mc_eff_vs_pt_TnP%d_PtW%d_cent_%i_to_%i_absy1p8_2p4",isTnP,isPtW,cLow,cHigh));}


  TFile *fAcc = new TFile(Form("./Eff_Acc/roots/acceptance_Prompt_GenOnly_wgt%d_210915.root",1),"read");
  TH1D* hAccPt[3];
  
  hAccPt[0] = (TH1D*) fAcc -> Get("hAccPt_2021_ally");
  hAccPt[1] = (TH1D*) fAcc -> Get("hAccPt_2021_midy");
  hAccPt[2] = (TH1D*) fAcc -> Get("hAccPt_2021_Fory");

  TFile *fFinal; TFile *fCErr;
  //fCErr = new TFile(Form("Macros/2021_09_14/roots/TnPSysUp_210928/2DFit_%s/CtauErr/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), "PR", fEffW, fAccW, isPtW, isTnP));

  double ctauLo;
  double ctauHi;
  double bfrac1;
  double bfrac2;
  double bfrac3;
  double ctauErrMin;
  double ctauErrMax;

  //if(ctauCut!=2){
  fCErr = new TFile(Form("Macros/2021_09_14/roots/2DFit_%s/CtauErr/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), "PR", 1, 1, 1, 1));
  if(ptLow==3&&ptHigh==6.5)fFinal = new TFile(Form("Macros/2021_09_14/roots/2DFit_%s/CtauRes_Sys/2DFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), "PR", 1, 1, 1, 1));
  else fFinal = new TFile(Form("Macros/2021_09_14/roots/2DFit_%s/CtauRes_Sys/2DFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), "PR", 1, 1, 1, 1));
  RooDataSet *dataw_Bkg = (RooDataSet*)fCErr->Get("dataw_Bkg");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataw_Bkg);
  if (ptLow==3 && ptHigh==6.5 && cLow==40){
    ctauErrMin = 0.02;//0.0143889
    ctauErrMax = 0.16;}
  else if (ptLow==7.5 && cLow==40){
    ctauErrMin = 0.015;//0.0143889
    ctauErrMax = 0.09;}
  else if (ptLow==9 && cLow==40){
    ctauErrMin = 0.015;//0.0143889
    ctauErrMax = 0.08;}//0.124872
  else if (ptLow==12 && cLow==40){
    ctauErrMin = 0.01;//0.0143889
    ctauErrMax = 0.05;}//0.124872
  else if (ptLow==6.5 && cLow==40 && cHigh==60){
    ctauErrMin = ws->var("ctau3DErr")->getMin();
    ctauErrMax = 0.112428;}//0.124872
  else {ctauErrMin = ws->var("ctau3DErr")->getMin();  ctauErrMax = ws->var("ctau3DErr")->getMax();}
  cout<<"CtauErr Min: "<<ctauErrMin<<", Max: "<<ctauErrMax<<endl;

  TH1D *Fraction1 = (TH1D*)fFinal->Get("Fraction1");
  TH1D *Fraction2 = (TH1D*)fFinal->Get("Fraction2");
  TH1D *Fraction3 = (TH1D*)fFinal->Get("Fraction3");

  ctauLo=Fraction1->GetBinLowEdge((double)Fraction1->FindFirstBinAbove(1e-3))+Fraction1->GetBinWidth((double)Fraction1->FindFirstBinAbove(1e-3));
  ctauHi=Fraction2->GetBinLowEdge((double)Fraction2->FindFirstBinAbove(1e-3))+Fraction2->GetBinWidth((double)Fraction2->FindFirstBinAbove(1e-3));
  bfrac1=Fraction1->GetBinContent((double)Fraction1->FindFirstBinAbove(1e-3));
  bfrac2=Fraction2->GetBinContent((double)Fraction2->FindFirstBinAbove(1e-3));
  bfrac3=Fraction3->GetBinContent((double)Fraction3->FindFirstBinAbove(1e-3));

  if(ctauCut==-1){
    cutName="ctauL"; cout<<"Ctau Cut: "<<cutName<<" < "<< ctauLo <<", "<< endl;}
  else if(ctauCut==0){
    cutName="ctauC"; cout<<"Ctau Cut: "<<cutName<<", "<< ctauLo <<"-"<< ctauHi<<endl;}
  else if(ctauCut==1){
    cutName="ctauR"; cout<<"Ctau Cut: "<<cutName<<" > "<< ctauHi <<", "<< endl;}
  else cutName="woCut";

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
  float ctau3DErr[nMaxDimu];
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
  TBranch *b_ctau3DErr;
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
  tree -> SetBranchAddress("ctau3DErr", ctau3DErr, &b_ctau3DErr);
  tree -> SetBranchAddress("qxa", qxa, &b_qxa);
  tree -> SetBranchAddress("qxb", qxb, &b_qxb);
  tree -> SetBranchAddress("qxc", qxc, &b_qxc);
  tree -> SetBranchAddress("qxdimu", qxdimu, &b_qxdimu);
  tree -> SetBranchAddress("qya", qya, &b_qya);
  tree -> SetBranchAddress("qyb", qyb, &b_qyb);
  tree -> SetBranchAddress("qyc", qyc, &b_qyc);
  tree -> SetBranchAddress("qydimu", qydimu, &b_qydimu);
  tree -> SetBranchAddress("weight", &weight, &b_weight);

  const int nMassBin = 9;

  float massBinDiff[nMassBin+1]={2.6, 2.8, 3.0, 3.06, 3.09, 3.12, 3.15, 3.2, 3.3, 3.5};
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
  TH1D* h_v3_1[nMass_div];
  TH1D* h_v3_2[nMass_div];
  TH1D* h_v3_3[nMass_div];
  TH1D* h_v3_4[nMass_div];
  TH1D* h_ljpsi[nMass_div];

  double Q_avg_low = -6500;
  double Q_avg_high = 27000;
  double Q_avg_low_dimu = -150;
  double Q_avg_high_dimu = 150;
  for(int imass=0; imass<nMass_div;imass++){
    h_v3_1[imass] = new TH1D(Form("h_v3_1_%d",imass),";#LTQ_{2}Q_{2A}^{*}#GT;counts",nQBin,Q_avg_low_dimu,Q_avg_high_dimu);
    h_v3_2[imass] = new TH1D(Form("h_v3_2_%d",imass),";#LTQ_{2A}Q_{2B}^{*}#GT;counts",nQBin,Q_avg_low,Q_avg_high);
    h_v3_3[imass] = new TH1D(Form("h_v3_3_%d",imass),";#LTQ_{2A}Q_{2C}^{*}#GT;counts",nQBin,Q_avg_low,Q_avg_high);
    h_v3_4[imass] = new TH1D(Form("h_v3_4_%d",imass),";#LTQ_{2B}Q_{2C}^{*}#GT;counts",nQBin,Q_avg_low,Q_avg_high);
    h_ljpsi[imass] = new TH1D(Form("h_ljBkg_%d",imass),";l_{J/#psi(mm)};Counts/(0.03 mm)",1050,-1.5,2);
  }

  TH1D* h_v3_num_q1 = new TH1D("h_v3_num_q1",";m_{#mu^{+}#mu^{-}};#LTQ_{2}Q_{2A}^{*}#GT" ,nMassBin,massBin);
  TH1D* h_v3_den_q2 = new TH1D("h_v3_num_q2",";m_{#mu^{+}#mu^{-}};#LTQ_{2A}Q_{2B}^{*}#GT",nMassBin,massBin);
  TH1D* h_v3_den_q3 = new TH1D("h_v3_num_q3",";m_{#mu^{+}#mu^{-}};#LTQ_{2A}Q_{2C}^{*}#GT",nMassBin,massBin);
  TH1D* h_v3_den_q4 = new TH1D("h_v3_num_q4",";m_{#mu^{+}#mu^{-}};#LTQ_{2B}Q_{2C}^{*}#GT",nMassBin,massBin);

  TH1D* h_mass = new TH1D("h_mass",";m_{#mu^{+}#mu^{-}};Counts/(0.025 GeV)",36,massLow,massHigh);
  TH1D* h_decayL = new TH1D("h_decayL",";l_{J/#psi(mm)};Counts/(0.1 mm)",100,-4,6);
  TH1D* h_pt = new TH1D("h_pt",";p_{T};Counts/",900,0,9);
  //TH1D* h_ljpsi_bkg = new TH1D("h_ljpsi_bkg",";l_{J/#psi(mm)};Counts/(0.03 mm)",1050,-1.5,2);

  RooRealVar* massVar = new RooRealVar("mass","mass variable",0,200,"GeV/c^{2}");
  RooRealVar* weightVar = new RooRealVar("weight", "weight var", 0, 10000,"");

  RooArgSet* argSet = new RooArgSet(*massVar, *weightVar);
  RooDataSet* dataSet = new RooDataSet("dataset","a dataset",*argSet); 

  const static int countMax = 1000000;
  vector<vector<double>> v3_1(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v3_2(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v3_3(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v3_4(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v3_1_raw(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v3_2_raw(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v3_3_raw(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> v3_4_raw(nMassBin,vector<double> (countMax,0));
  vector<vector<double>> weight_dimu(nMassBin,vector<double> (countMax,0));

  vector<double> v3_1_avg(nMassBin,0);
  vector<double> v3_2_avg(nMassBin,0);
  vector<double> v3_3_avg(nMassBin,0);
  vector<double> v3_4_avg(nMassBin,0);

  int dbcount=0;
  vector<unsigned int> count(nMassBin,0);
  vector<unsigned int> count_ss(nMassBin,0);

  vector<double> weight_s(nMassBin,0);

  double weight_acc = 1;
  double weight_eff = 1;

  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;

  int nDimu_all=0;
  int nDimuPass=0;
  int nDimu_one=0;
  int nDimu_more=0;
  //Begin Loop
  for(int i=0; i<nEvt; i++){
    tree->GetEntry(i);
    nDimuPass=0;
    if(cBin>cLow&&cBin<cHigh){ 
      if(fabs(vz)<15){
        for(int j=0; j<nDimu; j++){
          if(! ((double)pt[j]>3 && (double)pt[j]<50 && abs((double)y[j])<2.4 && IsAcceptanceQQ(pt1[j],eta1[j]) && IsAcceptanceQQ(pt2[j],eta2[j])) ) continue;
          nDimuPass++;
        }
        nDimu_all++;
        if(nDimuPass>1) {nDimu_more++; continue;}
        if(nDimuPass==1) nDimu_one++;
        // Fill Dimuon Loop
        for(int j=0; j<nDimu; j++){
          if(pt[j]>ptLow&&pt[j]<ptHigh&&recoQQsign[j]==0&&mass[j]>massLow&&mass[j]<massHigh&&abs(y[j])>yLow&&abs(y[j])<yHigh&&
              ctau3DErr[j]>ctauErrMin&&ctau3DErr[j]<ctauErrMax){
            if (ctauCut==-1 && !(ctau3D[j]<ctauLo)) continue;
            else if (ctauCut==0 && !(ctau3D[j]>ctauLo && ctau3D[j]<ctauHi)) continue;
            else if (ctauCut==1 && !(ctau3D[j]>ctauHi)) continue;
            weight_acc=1;
            weight_eff=1;
            if(isPtW){
              if ( abs((double)y[j])<1.6 ) {weight_acc = getAccWeight(hAccPt[1], pt[j]);}
                //cout <<"hAccPt1, "<<"Acc Weight : " << (double)1.0/weight_acc << " y : " << y[j] << " pT : " << pt[j] << endl;}
              else if ( abs((double)y[j])>1.6 && abs((double)y[j])<2.4 ) { weight_acc = getAccWeight(hAccPt[2], pt[j]);}
                //cout <<"hAccPt2, "<<"Acc Weight : " << (double)1.0/weight_acc << " y : " << y[j] << " pT : " << pt[j] << endl;}
            }
            if(fEffW){
              if ( abs((double)y[j])>=0.0 && abs((double)y[j])<1.2 ) { weight_eff = getEffWeight(hEffPt[0], pt[j]); }
              else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) { weight_eff = getEffWeight(hEffPt[1], pt[j]); }
              else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<1.8 ) { weight_eff = getEffWeight(hEffPt[2], pt[j]); }
              else if ( abs((double)y[j])>=1.8 && abs((double)y[j])<2.4 ) { weight_eff = getEffWeight(hEffPt[3], pt[j]); }
            }
            double weight_ = weight * weight_eff * weight_acc;
            for(int imbin=0; imbin<nMassBin; imbin++){
              if(mass[j]>=massBin[imbin] && mass[j]<massBin[imbin+1]){
                v3_1[imbin][count[imbin]] = (qxa[j]*qxdimu[j] + qya[j]*qydimu[j])*weight_;
                v3_2[imbin][count[imbin]] = (qxa[j]*qxb[j] + qya[j]*qyb[j])*weight_;
                v3_3[imbin][count[imbin]] = (qxa[j]*qxc[j] + qya[j]*qyc[j])*weight_;
                v3_4[imbin][count[imbin]] = (qxb[j]*qxc[j] + qyb[j]*qyc[j])*weight_;
              
                v3_1_raw[imbin][count[imbin]] = (qxa[j]*qxdimu[j] + qya[j]*qydimu[j]);
                v3_2_raw[imbin][count[imbin]] = (qxa[j]*qxb[j] + qya[j]*qyb[j]);
                v3_3_raw[imbin][count[imbin]] = (qxa[j]*qxc[j] + qya[j]*qyc[j]);
                v3_4_raw[imbin][count[imbin]] = (qxb[j]*qxc[j] + qyb[j]*qyc[j]);

                v3_1_avg[imbin] += v3_1[imbin][count[imbin]];
                v3_2_avg[imbin] += v3_2[imbin][count[imbin]];
                v3_3_avg[imbin] += v3_3[imbin][count[imbin]];
                v3_4_avg[imbin] += v3_4[imbin][count[imbin]];

                weight_dimu[imbin][count[imbin]] = weight_;

                count[imbin]++;
                weight_s[imbin] += weight_;
              }
            }
            if(mass[j]>=mass_low_SB1 && mass[j]<mass_high_SB1){
              h_v3_1[0]->Fill(qxa[j]*qxdimu[j] + qya[j]*qydimu[j], weight_);
              h_v3_2[0]->Fill(qxa[j]*qxb[j] + qya[j]*qyb[j], weight_);
              h_v3_3[0]->Fill(qxa[j]*qxc[j] + qya[j]*qyc[j], weight_);
              h_v3_4[0]->Fill(qxb[j]*qxc[j] + qyb[j]*qyc[j], weight_);
              h_ljpsi[0]->Fill(ctau3D[j],weight_);
            }
            else if(mass[j]>=mass_low_SB2 && mass[j]<mass_high_SB2){
              h_v3_1[1]->Fill(qxa[j]*qxdimu[j] + qya[j]*qydimu[j], weight_);
              h_v3_2[1]->Fill(qxa[j]*qxb[j] + qya[j]*qyb[j], weight_);
              h_v3_3[1]->Fill(qxa[j]*qxc[j] + qya[j]*qyc[j], weight_);
              h_v3_4[1]->Fill(qxb[j]*qxc[j] + qyb[j]*qyc[j], weight_);
              h_ljpsi[1]->Fill(ctau3D[j],weight_);
            }
            else if(mass[j]>=mass_low_SB && mass[j]<mass_high_SB){
              h_v3_1[2]->Fill(qxa[j]*qxdimu[j] + qya[j]*qydimu[j], weight_);
              h_v3_2[2]->Fill(qxa[j]*qxb[j] + qya[j]*qyb[j], weight_);
              h_v3_3[2]->Fill(qxa[j]*qxc[j] + qya[j]*qyc[j], weight_);
              h_v3_4[2]->Fill(qxb[j]*qxc[j] + qyb[j]*qyc[j], weight_);
              h_ljpsi[2]->Fill(ctau3D[j],weight_);
            }
            massVar->setVal((double)mass[j]);
            weightVar->setVal((double)weight_);
            dataSet->add(*argSet);
            h_mass->Fill(mass[j],weight_);
            if(mass[j]>=3.06&&mass[j]<=3.09)h_decayL->Fill(ctau3D[j]);
            if(mass[j]>=3.06&&mass[j]<=3.09)h_pt->Fill(pt[j]);
          }
        }
      }
    }
  }
  cout<<"How many Jpsi??: "<<float(h_mass->GetEntries())<<endl;

  TLine *lcut;

  int nMassFrameBin = 36;
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

  vector<double> v3_1_err(nMassBin,0);
  vector<double> v3_2_err(nMassBin,0);
  vector<double> v3_3_err(nMassBin,0);
  vector<double> v3_4_err(nMassBin,0);

  for(int ibin=0; ibin<nMassBin; ibin++){
    v3_1_avg[ibin] = v3_1_avg[ibin]/weight_s[ibin];
    v3_2_avg[ibin] = v3_2_avg[ibin]/weight_s[ibin];
    v3_3_avg[ibin] = v3_3_avg[ibin]/weight_s[ibin];
    v3_4_avg[ibin] = v3_4_avg[ibin]/weight_s[ibin];

    for(int icount=0; icount<count[ibin]; icount++){
      v3_1_err[ibin] += (v3_1_raw[ibin][icount]-v3_1_avg[ibin])*(v3_1_raw[ibin][icount]-v3_1_avg[ibin]) * weight_dimu[ibin][icount] * weight_dimu[ibin][icount];
      v3_2_err[ibin] += (v3_2_raw[ibin][icount]-v3_2_avg[ibin])*(v3_2_raw[ibin][icount]-v3_2_avg[ibin]) * weight_dimu[ibin][icount] * weight_dimu[ibin][icount];
      v3_3_err[ibin] += (v3_3_raw[ibin][icount]-v3_3_avg[ibin])*(v3_3_raw[ibin][icount]-v3_3_avg[ibin]) * weight_dimu[ibin][icount] * weight_dimu[ibin][icount];
      v3_4_err[ibin] += (v3_4_raw[ibin][icount]-v3_4_avg[ibin])*(v3_4_raw[ibin][icount]-v3_4_avg[ibin]) * weight_dimu[ibin][icount] * weight_dimu[ibin][icount];
    }

    v3_1_err[ibin] = TMath::Sqrt(v3_1_err[ibin])/weight_s[ibin];
    v3_2_err[ibin] = TMath::Sqrt(v3_2_err[ibin])/weight_s[ibin];
    v3_3_err[ibin] = TMath::Sqrt(v3_3_err[ibin])/weight_s[ibin];
    v3_4_err[ibin] = TMath::Sqrt(v3_4_err[ibin])/weight_s[ibin];

    h_v3_num_q1->SetBinContent(ibin+1,v3_1_avg[ibin]);
    h_v3_num_q1->SetBinError(ibin+1,v3_1_err[ibin]);
    h_v3_den_q2->SetBinContent(ibin+1,v3_2_avg[ibin]);
    h_v3_den_q2->SetBinError(ibin+1,v3_2_err[ibin]);
    h_v3_den_q3->SetBinContent(ibin+1,v3_3_avg[ibin]);
    h_v3_den_q3->SetBinError(ibin+1,v3_3_err[ibin]);
    h_v3_den_q4->SetBinContent(ibin+1,v3_4_avg[ibin]);
    h_v3_den_q4->SetBinError(ibin+1,v3_4_err[ibin]);

    cout << ibin << "th Bin : " << count[ibin] << ",  weight sum  : " << weight_s[ibin] << endl;
    cout << "v3_1_avg " << ibin << " : " << v3_1_avg[ibin] << endl;
    cout << "v3_2_avg " << ibin << " : " << v3_2_avg[ibin] << endl;
    cout << "v3_3_avg " << ibin << " : " << v3_3_avg[ibin] << endl;
    cout << "v3_4_avg " << ibin << " : " << v3_4_avg[ibin] << endl;

    cout << "h_v3_num_q1 " << ibin << ", val : " << h_v3_num_q1->GetBinContent(ibin+1) << " err : " << h_v3_num_q1->GetBinError(ibin+1) << endl;
    cout << "h_v3_den_q2 " << ibin << ", val : " << h_v3_den_q2->GetBinContent(ibin+1) << " err : " << h_v3_den_q2->GetBinError(ibin+1) << endl;
    cout << "h_v3_den_q3 " << ibin << ", val : " << h_v3_den_q3->GetBinContent(ibin+1) << " err : " << h_v3_den_q3->GetBinError(ibin+1) << endl;
    cout << "h_v3_den_q4 " << ibin << ", val : " << h_v3_den_q4->GetBinContent(ibin+1) << " err : " << h_v3_den_q4->GetBinError(ibin+1) << endl;
  }

  TH1D* h_v3_den_ = (TH1D*)h_v3_den_q2->Clone("h_v3_den_");
  h_v3_den_->Multiply(h_v3_den_q3);
  h_v3_den_->Divide(h_v3_den_q4);

  TH1D* h_v3_den = (TH1D*) h_v3_den_->Clone("h_v3_den"); //h_v3_den->Reset();
  GetHistSqrt(h_v3_den_,h_v3_den);

  TH1D* h_v3_final = (TH1D*) h_v3_num_q1 -> Clone("h_v3_SplusB");
  h_v3_final->Divide(h_v3_den);

  SetHistStyle(h_v3_final,0,0);
  SetGraphStyle2(g_mass,0,0);
  //g_mass->GetYaxis()->SetLimits(0,10000);
  //g_mass->GetYaxis()->SetRangeUser(0,);
  g_mass->GetYaxis()->SetTitle("Events/(0.025 GeV/c^{2})");
  g_mass->SetMarkerSize(1);
  g_mass->GetYaxis()->SetLabelSize(0.055);
  g_mass->GetXaxis()->SetLabelSize(0.04);
  g_mass->GetYaxis()->SetTitleSize(0.055);
  g_mass->GetYaxis()->SetTitleOffset(1.7);
  g_mass->GetXaxis()->SetTitleSize(0.065);

  h_v3_final->GetYaxis()->SetRangeUser(-0.05,0.26);
  h_v3_final->GetYaxis()->SetTitle("v_{2}^{S+B}");
  h_v3_final->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  h_v3_final->GetYaxis()->SetLabelSize(0.055);
  h_v3_final->GetXaxis()->SetLabelSize(0.055);
  h_v3_final->GetXaxis()->SetTitleSize(0.07);
  h_v3_final->GetYaxis()->SetTitleSize(0.07);
  h_v3_final->GetYaxis()->SetTitleOffset(1.2);
  h_v3_final->GetXaxis()->SetTitleOffset(1.0);
  h_v3_final->GetXaxis()->SetLabelOffset(0.011);
  TLine *l=new TLine(2.6, 0, 3.5, 0);
  l->SetLineColor(kBlack);

  float pos_x = 0.43;
  float pos_x_mass = 0.73;
  float pos_y = 0.7;
  float pos_y_diff = 0.071;
  int text_color = 1;
  float text_size = 16;
  TString perc = "%";

  //TGaxis::SetMaxDigits(2);
  TCanvas* c_mass_v3 = new TCanvas("c_mass_v3","",1300,0,590,750);
  TPad* pad1 = new TPad("pad1","pad1",0,0.5,1.0,1.0);
  TPad* pad2 = new TPad("pad2","pad2",0,0.0,1.0,0.5);
  c_mass_v3->cd();
  pad1->SetTicks(1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLeftMargin(0.19);
  pad1->SetTopMargin(0.08);
  pad1->cd();
  pad1->SetLogy();
  g_mass->GetXaxis()->SetTitleSize(0);
  g_mass->GetXaxis()->SetLabelSize(0);
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
  g_mass->GetXaxis()->SetLimits(massLow,massHigh);
  g_mass->GetXaxis()->SetRangeUser(massLow,massHigh);

  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  else if(ptLow ==6.5) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  else if(ptHigh==6.5) drawText(Form("%.f < p_{T}^{#mu#mu} < %.1f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  else if(ptLow!=0) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",(float)ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  //drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
  //    drawText(Form("#font[12]{l}_{J/#psi} < %.4f", ctauCut),pos_x_mass,pos_y-pos_y_diff*5,text_color,text_size);
  if(ctauCut==-1){
    drawText(Form("#font[12]{l}_{J/#psi} < %.5f", ctauLo),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);
    drawText(Form("NP J/#psi Frac. %.3f", bfrac1),pos_x_mass,pos_y-pos_y_diff*5,text_color,text_size);}
  else if(ctauCut==0){
    drawText(Form("%.5f < #font[12]{l}_{J/#psi} < %.5f", ctauLo, ctauHi),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);
    drawText(Form("NP J/#psi Frac. %.3f", bfrac3),pos_x_mass,pos_y-pos_y_diff*5,text_color,text_size);}
  else if(ctauCut==1){
    drawText(Form("#font[12]{l}_{J/#psi} > %.5f", ctauHi),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);
    drawText(Form("NP J/#psi Frac. %.3f", bfrac2),pos_x_mass,pos_y-pos_y_diff*5,text_color,text_size);}
  else drawText("#font[12]{l}_{J/#psi} Inclusive",pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.17);
  pad2->SetLeftMargin(0.19);
  pad2->SetTicks();
  pad2->cd();
  ////////////pad2//////////
  jumSun(massLow,0,massHigh,0,1,1);
  h_v3_final->Draw("P");
  l->Draw("same");

  //
//CMS_lumi_v2mass(pad1,iPeriod,iPos,1);  
  
CMS_lumi_v2mass(pad1,iPeriod,iPos);  
  pad1->Update();
  pad2->Update();
  c_mass_v3->cd();
  pad1->Draw();
  pad2->Draw();
  //gSystem->mkdir(Form("figs/TnPSysUp_210928/q_vector",sample.Data()),1);
  if(ctauCut==0) c_mass_v3->SaveAs(Form("figs/TnPSysUp_210928/v3mass_hist/%s/v3Mass_%s_ctau_%.5f_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauLo,ctauHi,cutName.Data()));
  else if(ctauCut==-1)c_mass_v3->SaveAs(Form("figs/TnPSysUp_210928/v3mass_hist/%s/v3Mass_%s_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauLo,cutName.Data()));
  else if(ctauCut== 1)c_mass_v3->SaveAs(Form("figs/TnPSysUp_210928/v3mass_hist/%s/v3Mass_%s_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauHi,cutName.Data()));
  else c_mass_v3->SaveAs(Form("figs/TnPSysUp_210928/v3mass_hist/%s/v3Mass_%s_Inc.pdf",DATE.Data(),kineLabel.Data()));

  TCanvas* c_pt = new TCanvas("c_pt","",0,500,600,600);
  //c_pt->SetLogy();
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  else if(ptLow!=0) drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  //drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
  h_pt->Draw("P");
  //lcut->Draw("same");
  c_pt->Update();
  if(ctauCut==0) c_pt->SaveAs(Form("figs/TnPSysUp_210928/pt/%s/pt_%s_ctau_%.5f_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauLo,ctauHi,cutName.Data()));
  else if(ctauCut==-1)c_pt->SaveAs(Form("figs/TnPSysUp_210928/pt/%s/pt_%s_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauLo,cutName.Data()));
  else if(ctauCut== 1)c_pt->SaveAs(Form("figs/TnPSysUp_210928/pt/%s/pt_%s_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauHi,cutName.Data()));
  else c_pt->SaveAs(Form("figs/TnPSysUp_210928/pt/%s/pt_%s_Inc.pdf",DATE.Data(),kineLabel.Data()));

  TCanvas* c_decayL = new TCanvas("c_decayL","",0,500,600,600);
  c_decayL->SetLogy();
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  else if(ptLow!=0) drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
  h_decayL->Draw("P");
  //lcut->Draw("same");
  c_decayL->Update();
  if(ctauCut==0) c_decayL->SaveAs(Form("figs/TnPSysUp_210928/decayL/%s/decayL_%s_ctau_%.5f_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauLo,ctauHi,cutName.Data()));
  else if(ctauCut==-1)c_decayL->SaveAs(Form("figs/TnPSysUp_210928/decayL/%s/decayL_%s_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauLo,cutName.Data()));
  else if(ctauCut== 1)c_decayL->SaveAs(Form("figs/TnPSysUp_210928/decayL/%s/decayL_%s_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauHi,cutName.Data()));
  else c_decayL->SaveAs(Form("figs/TnPSysUp_210928/decayL/%s/decayL_%s_Inc.pdf",DATE.Data(),kineLabel.Data()));

  TCanvas* c_mass = new TCanvas("c_mass","",600,600);
  c_mass->cd();
  h_mass->Draw("P");
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  else if(ptLow!=0) drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
  drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);
  
CMS_lumi_v2mass(c_mass,iPeriod,iPos);
  c_mass->Update();
  if(ctauCut==0)c_mass->SaveAs(Form("figs/TnPSysUp_210928/mass_dist/%s/MassDist_%s_ctau_%.5f_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauLo,ctauHi,cutName.Data()));
  else if(ctauCut==-1)c_mass->SaveAs(Form("figs/TnPSysUp_210928/mass_dist/%s/MassDist_%s_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauLo,cutName.Data()));
  else if(ctauCut==1)c_mass->SaveAs(Form("figs/TnPSysUp_210928/mass_dist/%s/MassDist_%s_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),ctauHi,cutName.Data()));
  else c_mass->SaveAs(Form("figs/TnPSysUp_210928/mass_dist/%s/MassDist_%s_Inc.pdf",DATE.Data(),kineLabel.Data()));

  TLegend *leg_v3_1 = new TLegend(0.38,0.64,0.77,0.9);
  TLegend *leg_v3_2 = new TLegend(0.38,0.64,0.77,0.9);
  TLegend *leg_v3_3 = new TLegend(0.38,0.64,0.77,0.9);
  TLegend *leg_v3_4 = new TLegend(0.38,0.64,0.77,0.9);
  SetLegendStyle(leg_v3_1);
  SetLegendStyle(leg_v3_2);
  SetLegendStyle(leg_v3_3);
  SetLegendStyle(leg_v3_4);

  double mean_v3_1[nMass_div];
  double mean_v3_2[nMass_div];
  double mean_v3_3[nMass_div];
  double mean_v3_4[nMass_div];

  double mean_err_v3_1[nMass_div];
  double mean_err_v3_2[nMass_div];
  double mean_err_v3_3[nMass_div];
  double mean_err_v3_4[nMass_div];

  for(int i=0;i<nMass_div;i++){
    SetHistStyle(h_v3_1[i],i,i);
    SetHistStyle(h_v3_2[i],i,i);
    SetHistStyle(h_v3_3[i],i,i);
    SetHistStyle(h_v3_4[i],i,i);

    scaleInt(h_v3_1[i]);
    scaleInt(h_v3_2[i]);
    scaleInt(h_v3_3[i]);
    scaleInt(h_v3_4[i]);

    mean_v3_1[i] = h_v3_1[i]->GetMean();
    mean_v3_2[i] = h_v3_2[i]->GetMean();
    mean_v3_3[i] = h_v3_3[i]->GetMean();
    mean_v3_4[i] = h_v3_4[i]->GetMean();

    mean_err_v3_1[i] = h_v3_1[i]->GetMeanError();
    mean_err_v3_2[i] = h_v3_2[i]->GetMeanError();
    mean_err_v3_3[i] = h_v3_3[i]->GetMeanError();
    mean_err_v3_4[i] = h_v3_4[i]->GetMeanError();

    leg_v3_1->AddEntry(h_v3_1[i],(fSB[i]+Form(" mean:%.2f",mean_v3_1[i])+Form(" err:%.2f",mean_err_v3_1[i])).Data(),"l");
    leg_v3_2->AddEntry(h_v3_2[i],(fSB[i]+Form(" mean:%.2f",mean_v3_2[i])+Form(" err:%.2f",mean_err_v3_2[i])).Data(),"l");
    leg_v3_3->AddEntry(h_v3_3[i],(fSB[i]+Form(" mean:%.2f",mean_v3_3[i])+Form(" err:%.2f",mean_err_v3_3[i])).Data(),"l");
    leg_v3_4->AddEntry(h_v3_4[i],(fSB[i]+Form(" mean:%.2f",mean_v3_4[i])+Form(" err:%.2f",mean_err_v3_4[i])).Data(),"l");
  }
  //cout<<"what?? : "<<h_v3_1[0]->GetEntries()<<endl;

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

  TCanvas *c_qq_1 = new TCanvas("c_qqa","",700,0,600,600);
  c_qq_1->cd();
  h_v3_1[2]->Draw("hist");
  h_v3_1[1]->Draw("hist same");
  h_v3_1[0]->Draw("hist same");
  leg_v3_1->Draw("same");
  leg_v3_1->SetFillStyle(0);
  if(ctauCut==0) c_qq_1->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qqa_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLo,ctauHi,cutName.Data()));
  else if(ctauCut==-1)c_qq_1->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qqa_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLo,cutName.Data()));
  else if(ctauCut==1)c_qq_1->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qqa_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauHi,cutName.Data()));
  else c_qq_1->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qqa_%s_Eff%d_Acc%d_PtW%d_TnP%d_Inc.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));

  TCanvas *c_qq_2 = new TCanvas("c_qaqb","",700,500,600,600);
  c_qq_2->cd();
  h_v3_2[2]->Draw("hist");
  h_v3_2[1]->Draw("hist same");
  h_v3_2[0]->Draw("hist same");
  leg_v3_2->Draw("same");
  if(ctauCut==0) c_qq_2->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qaqb_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLo,ctauHi,cutName.Data()));
  else if(ctauCut==-1)c_qq_2->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qaqb_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLo,cutName.Data()));
  else if(ctauCut==1)c_qq_2->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qaqb_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauHi,cutName.Data()));
  else c_qq_2->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qaqb_%s_Eff%d_Acc%d_PtW%d_TnP%d_Inc.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));

  TCanvas *c_qq_3 = new TCanvas("c_qaqc","",400,0,600,600);
  c_qq_3->cd();
  h_v3_3[2]->Draw("hist");
  h_v3_3[1]->Draw("hist same");
  h_v3_3[0]->Draw("hist same");
  leg_v3_3->Draw("same");
  if(ctauCut==0) c_qq_3->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qaqc_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLo,ctauHi,cutName.Data()));
  else if(ctauCut==-1)c_qq_3->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qaqc_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLo,cutName.Data()));
  else if(ctauCut== 1)c_qq_3->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qaqc_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauHi,cutName.Data()));
  else c_qq_3->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qaqc_%s_Eff%d_Acc%d_PtW%d_TnP%d_Inc.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));

  TCanvas *c_qq_4 = new TCanvas("c_qbqc","",400,500,600,600);
  c_qq_4->cd();
  h_v3_4[2]->Draw("hist");
  h_v3_4[1]->Draw("hist same");
  h_v3_4[0]->Draw("hist same");
  leg_v3_4->Draw("same");
  if(ctauCut==0) c_qq_4->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qbqc_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLo,ctauHi,cutName.Data()));
  else if(ctauCut==-1)c_qq_4->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qbqc_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLo,cutName.Data()));
  else if(ctauCut==1)c_qq_4->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qbqc_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauHi,cutName.Data()));
  else c_qq_4->SaveAs(Form("figs/TnPSysUp_210928/q_vector/%s/c_qbqc_%s_Eff%d_Acc%d_PtW%d_TnP%d_Inc.pdf",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP));

  TFile *wf; 
  if(ctauCut==0) wf = new TFile(Form("roots/TnPSysUp_210928/v3mass_hist/%s/Jpsi_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%.5f_%s.root",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLo,ctauHi,cutName.Data()),"recreate");
  else if(ctauCut==-1) wf = new TFile(Form("roots/TnPSysUp_210928/v3mass_hist/%s/Jpsi_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.root",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauLo,cutName.Data()),"recreate");
  else if(ctauCut==1) wf = new TFile(Form("roots/TnPSysUp_210928/v3mass_hist/%s/Jpsi_%s_Eff%d_Acc%d_PtW%d_TnP%d_ctau_%.5f_%s.root",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP,ctauHi,cutName.Data()),"recreate");
  else wf = new TFile(Form("roots/TnPSysUp_210928/v3mass_hist/%s/Jpsi_%s_Eff%d_Acc%d_PtW%d_TnP%d_Inc.root",DATE.Data(),kineLabel.Data(),fEffW,fAccW,isPtW,isTnP),"recreate");
  wf->cd();
  h_v3_final->Write();
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
  //cout<<"pt: "<<pt<<endl;
  double binN = h->FindBin(pt);
  double weight_ = 1./(h->GetBinContent(binN));
  return weight_;
}

double getEffWeight(TH1D *h, double pt){
  double binN = h->FindBin(pt);
  TF1 *eff1 = (TF1*)h->FindObject("f1");
  double eff = eff1->Eval(pt);
  double weight_ = 1./eff;
  return weight_;
}
