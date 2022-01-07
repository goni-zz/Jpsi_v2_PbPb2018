#include <iostream>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "JpsiUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBin.h"
#include "Style.h"
//#include "Style_jaebeom.h"
//#include "tnp_weight_lowptPbPb.h"
#include "tnp_weight_lowptPbPb_num_den_new.h"
#include <TAttMarker.h>

using namespace std;

void get_Eff_pt_psi_pbpb_dh_v3(
    float ptLow = 3.0, float ptHigh = 50.0,
    float yLow = 0.0, float yHigh = 2.4,
    float cLow = 20, float cHigh = 120, 
    //bool isTnP = false, bool isPtWeight = false, int state=1
    bool isTnP = true, bool isPtWeight = true, int isPtWgtUp = 0, // down 2, nominal 0
    int state=1
    ) {

  gStyle->SetOptStat(0);
  int kTrigSel = 12;	//jpsi=12,Upsilon=13

  float muPtCut = 0; //3.5, 1.8

  //jpsi
  float massLow = 0.0; 
  float massHigh = 10.0;
  //float massLow = 3.4;
  //float massHigh = 4.0;

  double min = 0;
  double max = ptHigh;
  const int numBins = 9; //50;//(max-min)/binwidth;  //31//


  //input files
  //PbPb
  TString inputMC1 = "../../OniatreeMC_Psi2S_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_pythia8.root";
  if(state==2) inputMC1 = "/work2/Oniatree/JPsi/OniatreeMC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8.root"; //PbPb_non prompt
  TChain* mytree = new TChain("myTree");
  mytree->Add(inputMC1.Data()); 
  //TFile *inf = new TFile("../OniatreeMC_Psi2S_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_pythia8.root","READ");
  //TTree *mytree = (TTree*)inf->Get("myTree");
  //if(state==1) mytree->Add(inputMC2.Data());
  cout<<"dmoon chk entries : "<<mytree->GetEntries()<<endl;
  //TString outFileName = "mc_eff_vs_pt_cent_prompt_pbpb_Jpsi.root"; 
  //if(state==2) outFileName = "mc_eff_vs_pt_cent_nprompt_pbpb_Jpsi.root"; 

  //pT reweighting function
  TFile *fPtW1 = new TFile("ratioDataMC_AA_Jpsi_DATA_All_y_211218.root","READ");
  TFile *fPtW2 = new TFile("ratioDataMC_AA_Jpsi_DATA_Forward_y_211218.root","read");
  if(state==2) TFile *fPtW = new TFile("ratioDataMC_AA_btojpsi_DATA_1s.root","read");
  TF1* fptw1 = (TF1*) fPtW1->Get("dataMC_Ratio1");
  TF1* fptw2 = (TF1*) fPtW2->Get("dataMC_Ratio1");

  // for systematics
  // parameters for all
  double p10 = 0.0, p11 = 0.0, p12 = 0.0, p13 = 0.0;
  // parameters for forward
  double p20 = 0.0, p21 = 0.0, p22 = 0.0, p23 = 0.0;
  // errors for all
  double e10 = 0.0, e11 = 0.0, e12 = 0.0, e13 = 0.0;
  // errors for forward
  double e20 = 0.0, e21 = 0.0, e22 = 0.0, e23 = 0.0;

  p10 = fptw1->GetParameter(0);
  p11 = fptw1->GetParameter(1);
  p12 = fptw1->GetParameter(2);
  p13 = fptw1->GetParameter(3);

  e10 = fptw1->GetParError(0);
  e11 = fptw1->GetParError(1);
  e12 = fptw1->GetParError(2);
  e13 = fptw1->GetParError(3);

  p20 = fptw2->GetParameter(0);
  p21 = fptw2->GetParameter(1);
  p22 = fptw2->GetParameter(2);
  p23 = fptw2->GetParameter(3);

  e20 = fptw2->GetParError(0);
  e21 = fptw2->GetParError(1);
  e22 = fptw2->GetParError(2);
  e23 = fptw2->GetParError(3);

  TF1 *fptw1sys = new TF1("fptw1sys","( [0]*x + [1]*x*x +[3]*x*x*x ) / (  (x-[2])*(x-[2])*(x-[2])  )",6.5, 50.0);
  TF1 *fptw2sys = new TF1("fptw2sys","( [0]*x + [1]*x*x +[3]*x*x*x ) / (  (x-[2])*(x-[2])*(x-[2])  )",3.0, 50.0);

  //cout<<"p10 : "<<p10<<", e10 : "<<e10<<", p20 : "<<p20<<", e20 : "<<e20<<endl;

  // sys for nominal
  if(isPtWgtUp == 0) {
    p10 = p10;
    p11 = p11;
    p12 = p12;
    p13 = p13;

    p20 = p20;
    p21 = p21;
    p22 = p22;
    p23 = p23;
  }

  // sys for up
  if(isPtWgtUp == 1) {
    p10 = p10 + e10;
    p11 = p11 + e11;
    p12 = p12 + e12;
    p13 = p13 + e13;

    p20 = p20 + e10;
    p21 = p21 + e11;
    p22 = p22 + e12;
    p23 = p23 + e13;
  }

  // sys for down
  if(isPtWgtUp == 2) {
    p10 = p10 - e10;
    p11 = p11 - e11;
    p12 = p12 - e12;
    p13 = p13 - e13;

    p20 = p20 - e10;
    p21 = p21 - e11;
    p22 = p22 - e12;
    p23 = p23 - e13;
  }

  //cout<<"p10 : "<<p10<<", e10 : "<<e10<<", p20 : "<<p20<<", e20 : "<<e20<<endl;

  fptw1sys->SetParameter(0, p10);
  fptw1sys->SetParameter(1, p11);
  fptw1sys->SetParameter(2, p12);
  fptw1sys->SetParameter(3, p13);

  fptw2sys->SetParameter(0, p20);
  fptw2sys->SetParameter(1, p21);
  fptw2sys->SetParameter(2, p22);
  fptw2sys->SetParameter(3, p23);


  double ptBin[numBins+1] = {0,3,4.5,6,7,8,9,10,15,50};
  float yBin[7] = {0,0.4,0.8,1.2,1.6,2.0,2.4};

  TH1D* hpt_reco_1 = new TH1D("hpt_reco_1","hpt_reco_1",numBins,ptBin);
  TH1D* hpt_reco_2 = new TH1D("hpt_reco_2","hpt_reco_2",numBins,ptBin);
  TH1D* hpt_reco_3 = new TH1D("hpt_reco_3","hpt_reco_3",numBins,ptBin);
  TH1D* hpt_reco_4 = new TH1D("hpt_reco_4","hpt_reco_4",numBins,ptBin);
  TH1D* hpt_reco_5 = new TH1D("hpt_reco_5","hpt_reco_5",numBins,ptBin);
  TH1D* hpt_reco_6 = new TH1D("hpt_reco_6","hpt_reco_5",numBins,ptBin);

  TH1D* hpt_gen_1 = new TH1D("hpt_gen_1","hpt_gen_1",numBins,ptBin);
  TH1D* hpt_gen_2 = new TH1D("hpt_gen_2","hpt_gen_2",numBins,ptBin);
  TH1D* hpt_gen_3 = new TH1D("hpt_gen_3","hpt_gen_3",numBins,ptBin);
  TH1D* hpt_gen_4 = new TH1D("hpt_gen_4","hpt_gen_4",numBins,ptBin);
  TH1D* hpt_gen_5 = new TH1D("hpt_gen_5","hpt_gen_5",numBins,ptBin);
  TH1D* hpt_gen_6 = new TH1D("hpt_gen_6","hpt_gen_6",numBins,ptBin);

  TF1 *f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])",3,50);
  f1->SetParameters(214,22,14);
  f1->SetParLimits(0,0,500);
  f1->SetParLimits(1,-1,500);
  f1->SetParLimits(2,0,500);

  TH1D* hy_gen = new TH1D("hy_gen","hy_gen",6,yBin);
  TH1D* hy_reco = new TH1D("hy_reco","hy_reco",6,yBin);
  hy_gen ->Sumw2();
  hy_reco ->Sumw2();

  hpt_reco_1->Sumw2();
  hpt_reco_2->Sumw2();
  hpt_reco_3->Sumw2();
  hpt_reco_4->Sumw2();
  hpt_reco_5->Sumw2();

  hpt_gen_1->Sumw2();
  hpt_gen_2->Sumw2();
  hpt_gen_3->Sumw2();
  hpt_gen_4->Sumw2();
  hpt_gen_5->Sumw2();

  const int maxBranchSize = 1000;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Int_t           Gen_QQ_size;
  Int_t           Gen_mu_size;
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_mu_4mom;
  ULong64_t       Gen_QQ_trig[maxBranchSize];   //[Gen_QQ_size]
  Float_t         Gen_QQ_VtxProb[maxBranchSize];   //[Gen_QQ_size]
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_mu_size;   //!
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Gen_mu_4mom;   //!
  TBranch        *b_Gen_QQ_trig;   //!
  TBranch        *b_Gen_QQ_VtxProb;   //!

  Gen_QQ_4mom = 0; Gen_mu_4mom = 0;
  mytree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  mytree->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
  mytree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  mytree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  //  muon id 
  Int_t           Gen_QQ_mupl_idx[maxBranchSize];
  Int_t           Gen_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Gen_QQ_mupl_idx;
  TBranch        *b_Gen_QQ_mumi_idx;
  mytree->SetBranchAddress("Gen_QQ_mupl_idx",Gen_QQ_mupl_idx,&b_Gen_QQ_mupl_idx);
  mytree->SetBranchAddress("Gen_QQ_mumi_idx",Gen_QQ_mumi_idx,&b_Gen_QQ_mumi_idx);

  Int_t           Gen_mu_charge[maxBranchSize];
  TBranch        *b_Gen_mu_charge;   //!
  mytree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);


  Int_t           Reco_QQ_size;
  Int_t           Reco_mu_size;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_mu_4mom;
  ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
  ULong64_t       Reco_mu_trig[maxBranchSize];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_mu_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_mu_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_mu_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  Reco_QQ_4mom = 0; Reco_mu_4mom = 0;
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

  //  muon id 
  Int_t           Reco_QQ_mupl_idx[maxBranchSize];
  Int_t           Reco_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Reco_QQ_mupl_idx;
  TBranch        *b_Reco_QQ_mumi_idx;
  mytree->SetBranchAddress("Reco_QQ_mupl_idx",Reco_QQ_mupl_idx,&b_Reco_QQ_mupl_idx);
  mytree->SetBranchAddress("Reco_QQ_mumi_idx",Reco_QQ_mumi_idx,&b_Reco_QQ_mumi_idx);


  Float_t         Reco_mu_dxy[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dxy;   //!
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  Float_t         Reco_mu_dz[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dz;   //!
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  Int_t           Reco_mu_nTrkWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Int_t           Reco_mu_nPixWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Int_t           Reco_QQ_sign[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);

  Int_t           Reco_mu_SelectionType[maxBranchSize];
  TBranch        *b_Reco_mu_SelectionType;
  mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

  Int_t Reco_mu_whichGen[maxBranchSize];
  TBranch *b_Reco_mu_whichGen;
  Float_t Gen_weight;
  TBranch *b_Gen_weight;
  mytree->SetBranchAddress("Reco_mu_whichGen",Reco_mu_whichGen, &b_Reco_mu_whichGen);
  mytree->SetBranchAddress("Gen_weight",&Gen_weight, &b_Gen_weight);


  TLorentzVector* JP_Gen= new TLorentzVector;
  TLorentzVector* mupl_Gen = new TLorentzVector;
  TLorentzVector* mumi_Gen = new TLorentzVector;

  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;

  double weight = 1;
  double tnp_weight = 1;
  double tnp_trig_weight = 1;
  double tnp_trig_weight_mupl = -1;
  double tnp_trig_weight_mumi = -1;
  double tnp_trig_weight_muplL2_num = -1;
  double tnp_trig_weight_muplL3_num = -1;
  double tnp_trig_weight_mumiL2_num = -1;
  double tnp_trig_weight_mumiL3_num = -1;
  double tnp_trig_weight_muplL2_den = -1;
  double tnp_trig_weight_muplL3_den = -1;
  double tnp_trig_weight_mumiL2_den = -1;
  double tnp_trig_weight_mumiL3_den = -1;
  double tnp_trig_weight_num = 1;
  double tnp_trig_weight_den = 1;
  double pt_weight = 1;

  double tnp_trig_dimu=-1;
  TH2D* hpt_tnp_trig = new TH2D("hpt_tnp_trig","hpt_tnp_trig",numBins,ptBin,40,0,2);

  int kL2filter = 16;	//jpsi=16,Upsilon=38
  int kL3filter = 17;	//jpsi=17,Upsilon=39

  int count =0;
  int counttnp =0;
  int nevt = mytree->GetEntries();
  //const int nevt = mytree->GetEntries();
  cout << "Total Events : " << nevt << endl;
  //nevt = 100;
  for(int iev=0; iev<nevt ; ++iev)
    //for(int iev=0; iev<300000 ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
    weight = findNcoll(Centrality) * Gen_weight;
    ///Gen_QQ_size = 1;
    //if(Gen_QQ_size > 0) cout<<"weight : "<<weight<<", Gen_QQ_size : "<<Gen_QQ_size<<", Reco_QQ_size : "<<Reco_QQ_size<<endl;

    for(int igen = 0; igen<Gen_QQ_size; igen++){
      JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(igen);
      mupl_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[igen]);
      mumi_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[igen]);

      Double_t Rapidity_g = fabs(JP_Gen->Rapidity());
      //cout<<"Rapidity_g : "<<Rapidity_g<<", pt : "<<JP_Gen->Pt()<<endl;
      if(! (JP_Gen->M()>massLow && JP_Gen->M()<massHigh) ) continue;

      if(! (JP_Gen->Pt() > 3 && JP_Gen->Pt()<50 && fabs(JP_Gen->Rapidity())<2.4 && IsAcceptanceQQ(mupl_Gen->Pt(),fabs(mupl_Gen->Eta()))&&IsAcceptanceQQ(mumi_Gen->Pt(),fabs(mumi_Gen->Eta()))) ) continue;

      if(! ( fabs(JP_Gen->Rapidity())<2.4 && fabs(mupl_Gen->Eta())<2.4 && fabs(mumi_Gen->Eta())<2.4) )  continue;
      if(Gen_mu_charge[Gen_QQ_mupl_idx[igen]]*Gen_mu_charge[Gen_QQ_mumi_idx[igen]]>0) continue;

      pt_weight = 1;
      //if(isPtWeight) pt_weight = fptw->Eval(JP_Gen->Pt());
      if(isPtWeight && fabs(JP_Gen->Rapidity()) < 2.4 && JP_Gen->Pt() > 6.5 ) pt_weight = fptw1sys->Eval(JP_Gen->Pt());
      if(isPtWeight && fabs(JP_Gen->Rapidity()) > 1.6 && JP_Gen->Pt() < 6.5) pt_weight = fptw2sys->Eval(JP_Gen->Pt());
      if(JP_Gen->Pt() > 6.5) hy_gen->Fill(Rapidity_g, weight*pt_weight);
      if(JP_Gen->Pt() > 6.5) hpt_gen_1->Fill(JP_Gen->Pt(),weight*pt_weight);

      if(! (Centrality > cLow && Centrality < cHigh)) continue;

      //cout<<"Rapidity_g : "<<Rapidity_g<<", pt : "<<JP_Gen->Pt()<<", weight : "<<weight<<", pt_weight : "<<pt_weight<<endl;
      if(Rapidity_g < 0.6 && JP_Gen->Pt() > 6.5) hpt_gen_2 -> Fill(JP_Gen->Pt(), weight*pt_weight);
      if(Rapidity_g > 0.0 && Rapidity_g < 1.2 && JP_Gen->Pt() > 6.5) hpt_gen_3 -> Fill(JP_Gen->Pt(), weight*pt_weight);
      if(Rapidity_g > 1.2 && Rapidity_g < 1.6 && JP_Gen->Pt() > 6.5) hpt_gen_4 -> Fill(JP_Gen->Pt(), weight*pt_weight);
      if(Rapidity_g > 1.6 && Rapidity_g < 1.8 && JP_Gen->Pt() > 3.0) hpt_gen_5 -> Fill(JP_Gen->Pt(), weight*pt_weight);
      if(Rapidity_g > 1.8 && Rapidity_g < 2.4 && JP_Gen->Pt() > 3.0) hpt_gen_6 -> Fill(JP_Gen->Pt(), weight*pt_weight);
    }


    bool HLTPass = false;
    if((HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) HLTPass=true;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
      //cout<<"Reco QQ pt : "<<JP_Reco->Pt()<<endl;


      bool HLTFilterPass=false;
      if( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) HLTFilterPass=true;

      if(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
      if(Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;

      bool passMuonTypePl = true;
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,1)));
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,3)));

      bool passMuonTypeMi = true;
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,1)));
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,3)));

      bool muplSoft = (  //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]])<20.) &&
          passMuonTypePl        //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
          ) ; 

      bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20.)  && 
          passMuonTypeMi       //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
          ) ; 

      if ( !(muplSoft && mumiSoft) ) continue;   
      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) continue;
      if(Reco_QQ_sign[irqq]!=0) continue;  

      Double_t Rapidity = fabs(JP_Reco->Rapidity());

      if(! (JP_Reco->Pt()>3&&JP_Reco->Pt()<50&&fabs(JP_Reco->Rapidity())<2.4&&IsAcceptanceQQ(mupl_Reco->Pt(),fabs(mupl_Reco->Eta()))&&IsAcceptanceQQ(mumi_Reco->Pt(),fabs(mumi_Reco->Eta()))) ) continue;

      if(!( fabs(mupl_Reco->Eta())<2.4 && fabs(mumi_Reco->Eta())<2.4 && fabs(JP_Reco->Rapidity())<2.4 && JP_Reco->M()>massLow && JP_Reco->M()<massHigh)) continue;


      if(HLTPass==true && HLTFilterPass==true) count++;
      if(isTnP){
        tnp_weight = 1;
        tnp_trig_weight = 1;
        tnp_trig_weight_mupl = -1;
        tnp_trig_weight_mumi = -1;
        tnp_trig_weight_muplL2_num = -1;
        tnp_trig_weight_muplL3_num = -1;
        tnp_trig_weight_mumiL2_num = -1;
        tnp_trig_weight_mumiL3_num = -1;
        tnp_trig_weight_muplL2_den = -1;
        tnp_trig_weight_muplL3_den = -1;
        tnp_trig_weight_mumiL2_den = -1;
        tnp_trig_weight_mumiL3_den = -1;
        tnp_trig_weight_num = 1;
        tnp_trig_weight_den = 1;


        tnp_weight = tnp_weight * tnp_weight_muid_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0) * tnp_weight_muid_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0); //mu id
        tnp_weight = tnp_weight * tnp_weight_trk_pbpb(mupl_Reco->Eta(), 0) * tnp_weight_trk_pbpb(mumi_Reco->Eta(), 0); //inner tracker

        //Trigger part
        if(!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ){
          //         cout << "irqq : " << irqq << " - iev : " << iev << endl;
          //         cout << "TnP ERROR !!!! ::: No matched L2 filter1 " << endl;
          continue;
        }
        bool mupl_L2Filter = ( (Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ? true : false ;
        bool mupl_L3Filter = ( (Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)) ) ? true : false ;
        bool mumi_L2Filter = ( (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ? true : false ;
        bool mumi_L3Filter = ( (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)) ) ? true : false ;
        if(mupl_L2Filter == false || mumi_L2Filter == false){ cout << "TnP ERROR !!!! ::: No matched L2 filter2 " << endl; cout << endl;} 

        bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter) ? true : false;
        bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter) ? true : false;
        bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter) ? true : false;
        bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter) ? true : false;
        bool SelDone = false;

        //if( mupl_isL2 && mumi_isL3){
        //        tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
        //        tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
        //        SelDone = true;
        //}
        //else if( mupl_isL3 && mumi_isL2){
        //        tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
        //        tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
        //        SelDone = true;
        //}
        //else if( mupl_isL3 && mumi_isL3){
        //        int t[2] = {-1,1}; // mupl, mumi
        //        int l = rand() % (2); 
        //        //pick up what will be L2
        //        if(t[l]==-1){
        // 	       tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
        // 	       tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
        //        }
        //        else if(t[l]==1){
        // 	       tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
        // 	       tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
        //        }
        //        else {cout << "ERROR :: No random selection done !!!!" << endl; continue;}
        //        SelDone = true;
        //}    

        if( mupl_isL2 && mumi_isL3){
          tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
          tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);
          SelDone = true;
          tnp_trig_weight = tnp_trig_weight_mupl * tnp_trig_weight_mumi;

        }
        else if( mupl_isL3 && mumi_isL2){
          tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
          tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
          SelDone = true;
          tnp_trig_weight = tnp_trig_weight_mupl * tnp_trig_weight_mumi;
        }
        else if( mupl_isL3 && mumi_isL3){

          tnp_trig_weight_muplL2_num = tnp_weight_trg_pbpb_num(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
          tnp_trig_weight_muplL3_num = tnp_weight_trg_pbpb_num(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
          tnp_trig_weight_mumiL2_num = tnp_weight_trg_pbpb_num(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
          tnp_trig_weight_mumiL3_num = tnp_weight_trg_pbpb_num(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);

          tnp_trig_weight_muplL2_den = tnp_weight_trg_pbpb_den(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
          tnp_trig_weight_muplL3_den = tnp_weight_trg_pbpb_den(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
          tnp_trig_weight_mumiL2_den = tnp_weight_trg_pbpb_den(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
          tnp_trig_weight_mumiL3_den = tnp_weight_trg_pbpb_den(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);

          tnp_trig_weight_num = tnp_trig_weight_muplL2_num * tnp_trig_weight_mumiL3_num + tnp_trig_weight_mumiL2_num * tnp_trig_weight_muplL3_num - tnp_trig_weight_muplL3_num * tnp_trig_weight_mumiL3_num;
          tnp_trig_weight_den = tnp_trig_weight_muplL2_den * tnp_trig_weight_mumiL3_den + tnp_trig_weight_mumiL2_den * tnp_trig_weight_muplL3_den - tnp_trig_weight_muplL3_den * tnp_trig_weight_mumiL3_den;
          tnp_trig_weight = tnp_trig_weight_num/tnp_trig_weight_den;
        }

        tnp_weight = tnp_weight * tnp_trig_weight;



        //if(SelDone == false){cout << "ERROR :: No muon filter combination selected !!!!" << endl; continue;}
        //if((tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1)){cout << "ERROR :: No trigger muon tnp scale factors selected !!!!" << endl; continue;}
        if(HLTPass==true && HLTFilterPass==true){
          counttnp++;
          tnp_trig_dimu = tnp_trig_weight;
          hpt_tnp_trig->Fill(JP_Reco->Pt(),tnp_trig_dimu);
        }
      }

      pt_weight = 1;
      //if(isPtWeight) pt_weight = fptw->Eval(JP_Reco->Pt());
      if(isPtWeight && fabs(JP_Reco->Rapidity()) < 2.4 && JP_Reco->Pt() > 6.5) pt_weight = fptw1sys->Eval(JP_Reco->Pt());
      if(isPtWeight && fabs(JP_Reco->Rapidity()) > 1.6 && JP_Reco->Pt() < 6.5) pt_weight = fptw2sys->Eval(JP_Reco->Pt());

      if(HLTPass==true && HLTFilterPass==true){
        if(JP_Reco->Pt() > 6.5 ) hy_reco->Fill(Rapidity, weight* tnp_weight* pt_weight);
        if(JP_Reco->Pt() > 6.5 ) hpt_reco_1->Fill(JP_Reco->Pt(),weight* tnp_weight* pt_weight);

        if(! (Centrality > cLow && Centrality < cHigh)) continue;

        if(Rapidity < 0.6 && JP_Reco->Pt() > 6.5) hpt_reco_2 -> Fill(JP_Reco->Pt(), weight* tnp_weight* pt_weight);
        if(Rapidity > 0.0 && Rapidity < 1.2 && JP_Reco->Pt() > 6.5) hpt_reco_3 -> Fill(JP_Reco->Pt(), weight* tnp_weight* pt_weight);
        if(Rapidity > 1.2 && Rapidity < 1.6 && JP_Reco->Pt() > 6.5) hpt_reco_4 -> Fill(JP_Reco->Pt(), weight* tnp_weight* pt_weight);
        if(Rapidity > 1.6 && Rapidity < 1.8 && JP_Reco->Pt() > 3.0) hpt_reco_5 -> Fill(JP_Reco->Pt(), weight* tnp_weight* pt_weight);
        if(Rapidity > 1.8 && Rapidity < 2.4 && JP_Reco->Pt() > 3.0) hpt_reco_6 -> Fill(JP_Reco->Pt(), weight* tnp_weight* pt_weight);
      }
    }
  }

  cout << "count " << count << endl;
  cout << "counttnp " << counttnp << endl;


  //Divide
  TH1D* hpt_eff_1;
  TH1D* hpt_eff_2;
  TH1D* hpt_eff_3;
  TH1D* hpt_eff_4;
  TH1D* hpt_eff_5;
  TH1D* hpt_eff_6;


  hpt_eff_1 = (TH1D*)hpt_reco_1->Clone("hpt_eff_1");
  hpt_eff_2 = (TH1D*)hpt_reco_2->Clone("hpt_eff_2");
  hpt_eff_3 = (TH1D*)hpt_reco_3->Clone("hpt_eff_3");
  hpt_eff_4 = (TH1D*)hpt_reco_4->Clone("hpt_eff_4");
  hpt_eff_5 = (TH1D*)hpt_reco_5->Clone("hpt_eff_5");
  hpt_eff_6 = (TH1D*)hpt_reco_6->Clone("hpt_eff_6");


  hpt_eff_1->Divide(hpt_eff_1, hpt_gen_1, 1, 1, "B");
  hpt_eff_2->Divide(hpt_eff_2, hpt_gen_2, 1, 1, "B");
  hpt_eff_3->Divide(hpt_eff_3, hpt_gen_3, 1, 1, "B");
  hpt_eff_4->Divide(hpt_eff_4, hpt_gen_4, 1, 1, "B");
  hpt_eff_5->Divide(hpt_eff_5, hpt_gen_5, 1, 1, "B");
  hpt_eff_6->Divide(hpt_eff_6, hpt_gen_6, 1, 1, "B");


  TH1D* hy_eff;
  hy_eff = (TH1D*)hy_reco->Clone("hy_eff");
  hy_eff->Divide(hy_eff, hy_gen, 1, 1, "B");

  hpt_eff_1->SetTitle("Eff: Rapidity 0.0-2.4");
  hpt_eff_2->SetTitle("Eff: Rapidity 0.0-0.6");
  hpt_eff_3->SetTitle("Eff: Rapidity 0.0-1.2");
  hpt_eff_4->SetTitle("Eff: Rapidity 1.2-1.6");
  hpt_eff_5->SetTitle("Eff: Rapidity 1.6-1.8");
  hpt_eff_6->SetTitle("Eff: Rapidity 1.8-2.4");

  f1->SetLineColor(kBlack);
  f1->SetLineWidth(2);
  hpt_eff_1->Fit(f1);
  f1->SetLineColor(kRed+1);
  f1->SetLineWidth(2);
  hpt_eff_2->Fit(f1);
  f1->SetLineColor(kOrange+1);
  f1->SetLineWidth(2);
  hpt_eff_3->Fit(f1);
  f1->SetLineColor(kGreen+1);
  f1->SetLineWidth(2);
  hpt_eff_4->Fit(f1);
  f1->SetLineColor(kBlue+1);
  f1->SetLineWidth(2);
  hpt_eff_5->Fit(f1);
  f1->SetLineColor(kViolet+1);
  f1->SetLineWidth(2);
  hpt_eff_6->Fit(f1);


  //draw same
  //gROOT->Macro("~/rootlogon.C");
  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.03);
  auto legend = new TLegend(0.6,0.84);

  gStyle->SetOptFit(0);
  TCanvas * cpt_eff = new TCanvas("cpt_eff","cpt_eff",0,0,900,800);
  cpt_eff->cd();
  hpt_eff_1->SetMarkerStyle(24);
  hpt_eff_1->SetMarkerColor(1);
  hpt_eff_1->SetLineColor(1);
  hpt_eff_1->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hpt_eff_1->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_1->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_1->Draw("E");
  //eeeeeee
  //hpt_eff_2->SetMarkerStyle(24);
  //hpt_eff_2->SetMarkerColor(1);
  //hpt_eff_2->SetLineColor(1);
  //hpt_eff_2->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  //hpt_eff_2->GetYaxis()->SetTitle("Efficiency");
  //hpt_eff_2->GetYaxis()->SetRangeUser(0.,1.3);
  //hpt_eff_2->Draw("E");
  //eeeeeee
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");
  //drawsame
  hpt_eff_2->SetMarkerStyle(25);
  hpt_eff_2->SetMarkerColor(kRed+1);
  hpt_eff_2->SetLineColor(kRed+1);
  hpt_eff_2->Draw("same");
  hpt_eff_3->SetMarkerStyle(30);
  hpt_eff_3->SetMarkerColor(kOrange+1);
  hpt_eff_3->SetLineColor(kOrange+1);
  hpt_eff_3->Draw("same");
  hpt_eff_4->SetMarkerStyle(26);
  hpt_eff_4->SetMarkerColor(kGreen+1);
  hpt_eff_4->SetLineColor(kGreen+1);
  hpt_eff_4->Draw("same");
  hpt_eff_5->SetMarkerStyle(27);
  hpt_eff_5->SetMarkerColor(kBlue+1);
  hpt_eff_5->SetLineColor(kBlue+1);
  hpt_eff_5->Draw("same");
  hpt_eff_6->SetMarkerStyle(28);
  hpt_eff_6->SetMarkerColor(kViolet+1);
  hpt_eff_6->SetLineColor(kViolet+1);
  hpt_eff_6->Draw("same");
  legend->AddEntry("hpt_eff_1",Form("|y|: 0-2.4, %0.0f-%0.0f %%",cLow/2,cHigh/2),"lep");
  legend->AddEntry("hpt_eff_2",Form("|y|: 0-1.2, %0.0f-%0.0f %%",cLow/2,cHigh/2),"lep");
  legend->AddEntry("hpt_eff_3",Form("|y|: 1.2-1.6, %0.0f-%0.0f %%",cLow/2,cHigh/2),"lep");
  legend->AddEntry("hpt_eff_4",Form("|y|: 1.6-1.8, %0.0f-%0.0f %%",cLow/2,cHigh/2),"lep");
  legend->AddEntry("hpt_eff_5",Form("|y|: 1.8-2.4, %0.0f-%0.0f %%",cLow/2,cHigh/2),"lep");
  legend->AddEntry("hpt_eff_6",Form("|y|: 1.6-2.4, %0.0f-%0.0f %%",cLow/2,cHigh/2),"lep");
  legend->SetBorderSize( 0);
  legend->Draw("E");
  cpt_eff->SaveAs("Eff_pt_noweight2.png");

  TCanvas * cy_eff = new TCanvas("cy_eff","cy_eff",0,0,900,800);
  cy_eff->cd();
  hy_eff->SetMarkerStyle(24);
  hy_eff->SetMarkerColor(1);
  hy_eff->SetLineColor(1);
  hy_eff->GetXaxis()->SetTitle("|y|");
  hy_eff->GetYaxis()->SetTitle("Efficiency");
  hy_eff->GetYaxis()->SetRangeUser(0.,1.2);
  hy_eff->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");
  cy_eff->SaveAs("Eff_absy_noweight2.png");

  //Save efficiency files for later use.

  hpt_eff_1 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0_2p4",isTnP, isPtWeight, cLow, cHigh));
  hpt_eff_2 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0_0p6  ",isTnP, isPtWeight, cLow, cHigh));
  hpt_eff_3 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy0p0_1p2",isTnP, isPtWeight, cLow, cHigh));
  hpt_eff_4 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy1p2_1p6",isTnP, isPtWeight, cLow, cHigh));
  hpt_eff_5 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy1p6_1p8",isTnP, isPtWeight, cLow, cHigh));
  hpt_eff_6 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_%0.0f_to_%0.0f_absy1p8_2p4",isTnP, isPtWeight, cLow, cHigh));
  hy_eff ->SetName(Form("mc_eff_vs_rap_TnP%d_PtW%d",isTnP, isPtWeight));

  //TString outFileName = Form("mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_prompt_pbpb_Jpsi_PtW%d_tnp%d_drawsame1.root",cLow,cHigh,isPtWeight,isTnP);
  TString outFileName = Form("mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_prompt_pbpb_psi2s_PtW%d_tnp%d_SysUp%d_20211218.root",cLow,cHigh,isPtWeight,isTnP,isPtWgtUp);
  if(state==2) outFileName = Form("mc_eff_vs_pt_cent_%0.0f_to_%0.0f_rap_nprompt_pbpb_psi2s_PtW%d_tnp%d_drawsame.root",cLow,cHigh,isPtWeight,isTnP);
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hpt_eff_1->Write();
  hpt_eff_2->Write();
  hpt_eff_3->Write();
  hpt_eff_4->Write();
  hpt_eff_5->Write();
  hpt_eff_6->Write();

  hy_eff->Write();
  cpt_eff->Write();
  if(isTnP) hpt_tnp_trig->Write();

  hpt_reco_1->Write();
  hpt_reco_2->Write();
  hpt_reco_3->Write();
  hpt_reco_4->Write();
  hpt_reco_5->Write();
  hpt_reco_6->Write();

  hpt_gen_1->Write();
  hpt_gen_2->Write();
  hpt_gen_3->Write();
  hpt_gen_4->Write();
  hpt_gen_5->Write();
  hpt_gen_6->Write();

  fptw1sys->Write();
  fptw2sys->Write();

  hy_gen->Write();
  hy_reco->Write();

  outFile->Close();

}
