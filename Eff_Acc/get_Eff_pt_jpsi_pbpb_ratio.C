#include <iostream>

#include <TLorentzVector.h>
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../HiEvtPlaneList.h"
#include "../cutsAndBin.h"
#include "../Style_jaebeom.h"
#include "tnp_weight_lowptPbPb_num_den.h"
//#include "tnp_weight_lowptPbPb.h"
#include <TAttMarker.h>

using namespace std;

void get_Eff_pt_jpsi_pbpb_ratio( 
  float ptLow = 0.0, float ptHigh = 50.0,
  float yLow = 0.0, float yHigh = 2.4,
  int cLow = 0, int cHigh = 180, bool isTnP = true, bool isPtWeight = true, int state=2
  ) {

  gStyle->SetOptStat(0);
  int kTrigSel = 12;	//jpsi=12,Upsilon=13

  float muPtCut = 0; //3.5, 1.8

  //jpsi
  float massLow = 2.6;
  float massHigh = 3.5;

  double min = 0;
  double max = ptHigh;
  //const int numBins = 6; //50;//(max-min)/binwidth;  //31//
  const int numBins = 9; //50;//(max-min)/binwidth;  //31//

  //input files
  //PbPb
  TString inputMC = "/work2/Oniatree/JPsi/OniatreeMC_JPsi_pThat2_TunedCP5_HydjetDrumMB_5p02TeV_Pythia8_part*.root";	//PbPb_prompt
  if(state==2) inputMC = "/work2/Oniatree/JPsi/OniatreeMC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8.root";	//PbPb_non prompt
  TChain* mytree = new TChain("myTree"); 
  mytree->Add(inputMC.Data());

  //TString outFileName = "mc_eff_vs_pt_cent_rap_prompt_pbpb_Jpsi_PtW.root"; 
  //if(state==2) outFileName = "mc_eff_vs_pt_cent_rap_nprompt_pbpb_Jpsi_PtW.root"; 

  //pT reweighting function
  TFile *fPtW = new TFile("WeightedFcN_fit/ratioDataMC_AA_DATA_1s.root","read");
  if(state==2) TFile *fPtW = new TFile("WeightedFcN_fit/ratioDataMC_AA_btojpsi_DATA_1s.root","read");
  TF1* fptw = (TF1*) fPtW->Get("dataMC_Ratio1");

  
  //double ptBin[numBins+1] = {0,3,4.5,6.5,8,12,50};
  double ptBin[numBins+1] = {0,3,4.5,6,7,8,9,10,15,50};

  //cent0-20
  TH1D* hpt_reco_1 = new TH1D("hpt_reco_1","hpt_reco_1",numBins,ptBin);
  TH1D* hpt_reco_2 = new TH1D("hpt_reco_2","hpt_reco_2",numBins,ptBin);
  TH1D* hpt_reco_3 = new TH1D("hpt_reco_3","hpt_reco_3",numBins,ptBin);
  TH1D* hpt_reco_4 = new TH1D("hpt_reco_4","hpt_reco_4",numBins,ptBin);

  //cent20-100
  TH1D* hpt_reco_5 = new TH1D("hpt_reco_5","hpt_reco_5",numBins,ptBin);
  TH1D* hpt_reco_6 = new TH1D("hpt_reco_6","hpt_reco_6",numBins,ptBin);
  TH1D* hpt_reco_7 = new TH1D("hpt_reco_7","hpt_reco_7",numBins,ptBin);
  TH1D* hpt_reco_8 = new TH1D("hpt_reco_8","hpt_reco_8",numBins,ptBin);

  //cent100-180
  TH1D* hpt_reco_9 = new TH1D("hpt_reco_9","hpt_reco_9",numBins,ptBin);
  TH1D* hpt_reco_10 = new TH1D("hpt_reco_10","hpt_reco_10",numBins,ptBin);
  TH1D* hpt_reco_11 = new TH1D("hpt_reco_11","hpt_reco_11",numBins,ptBin);
  TH1D* hpt_reco_12 = new TH1D("hpt_reco_12","hpt_reco_12",numBins,ptBin);

  //cent0-100
  TH1D* hpt_reco_13 = new TH1D("hpt_reco_13","hpt_reco_13",numBins,ptBin);
  TH1D* hpt_reco_14 = new TH1D("hpt_reco_14","hpt_reco_14",numBins,ptBin);
  TH1D* hpt_reco_15 = new TH1D("hpt_reco_15","hpt_reco_15",numBins,ptBin);
  TH1D* hpt_reco_16 = new TH1D("hpt_reco_16","hpt_reco_16",numBins,ptBin);
  TH1D* hpt_reco_17 = new TH1D("hpt_reco_17","hpt_reco_17",numBins,ptBin);

  //cent0-20
  TH1D* hpt_gen_1 = new TH1D("hpt_gen_1","hpt_gen_1",numBins,ptBin);
  TH1D* hpt_gen_2 = new TH1D("hpt_gen_2","hpt_gen_2",numBins,ptBin);
  TH1D* hpt_gen_3 = new TH1D("hpt_gen_3","hpt_gen_3",numBins,ptBin);
  TH1D* hpt_gen_4 = new TH1D("hpt_gen_4","hpt_gen_4",numBins,ptBin);

  //cent20-100
  TH1D* hpt_gen_5 = new TH1D("hpt_gen_5","hpt_gen_5",numBins,ptBin);
  TH1D* hpt_gen_6 = new TH1D("hpt_gen_6","hpt_gen_6",numBins,ptBin);
  TH1D* hpt_gen_7 = new TH1D("hpt_gen_7","hpt_gen_7",numBins,ptBin);
  TH1D* hpt_gen_8 = new TH1D("hpt_gen_8","hpt_gen_8",numBins,ptBin);
                                                   
  //cent100-180                                    
  TH1D* hpt_gen_9 = new TH1D("hpt_gen_9","hpt_gen_9",numBins,ptBin);
  TH1D* hpt_gen_10 = new TH1D("hpt_gen_10","hpt_gen_10",numBins,ptBin);
  TH1D* hpt_gen_11 = new TH1D("hpt_gen_11","hpt_gen_11",numBins,ptBin);
  TH1D* hpt_gen_12 = new TH1D("hpt_gen_12","hpt_gen_12",numBins,ptBin);

  //cent0-100
  TH1D* hpt_gen_13 = new TH1D("hpt_gen_13","hpt_gen_13",numBins,ptBin);
  TH1D* hpt_gen_14 = new TH1D("hpt_gen_14","hpt_gen_14",numBins,ptBin);
  TH1D* hpt_gen_15 = new TH1D("hpt_gen_15","hpt_gen_15",numBins,ptBin);
  TH1D* hpt_gen_16 = new TH1D("hpt_gen_16","hpt_gen_16",numBins,ptBin);
  TH1D* hpt_gen_17 = new TH1D("hpt_gen_17","hpt_gen_17",numBins,ptBin);


  TF1 *f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])",3,50);
  f1->SetParameters(214,22,14);
  f1->SetParLimits(0,0,500);
  f1->SetParLimits(1,-1,500);
  f1->SetParLimits(2,0,500);


  hpt_reco_1->Sumw2();
  hpt_reco_2->Sumw2();
  hpt_reco_3->Sumw2();
  hpt_reco_4->Sumw2();

  hpt_reco_5->Sumw2();
  hpt_reco_6->Sumw2();
  hpt_reco_7->Sumw2();
  hpt_reco_8->Sumw2();
            
  hpt_reco_9 ->Sumw2();
  hpt_reco_10->Sumw2();
  hpt_reco_11->Sumw2();
  hpt_reco_12->Sumw2();

  hpt_reco_13->Sumw2();
  hpt_reco_14->Sumw2();
  hpt_reco_15->Sumw2();
  hpt_reco_16->Sumw2();
  hpt_reco_17->Sumw2();
           
  hpt_gen_1->Sumw2();
  hpt_gen_2->Sumw2();
  hpt_gen_3->Sumw2();
  hpt_gen_4->Sumw2();

  hpt_gen_5->Sumw2();
  hpt_gen_6->Sumw2();
  hpt_gen_7->Sumw2();
  hpt_gen_8->Sumw2();
           
  hpt_gen_9->Sumw2();
  hpt_gen_10->Sumw2();
  hpt_gen_11->Sumw2();
  hpt_gen_12->Sumw2();

  hpt_gen_13->Sumw2();
  hpt_gen_14->Sumw2();
  hpt_gen_15->Sumw2();
  hpt_gen_16->Sumw2();
  hpt_gen_17->Sumw2();
           


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
  int countPtW =0;
  const int nevt = mytree->GetEntries();
  cout << "Total Events : " << nevt << endl;
  for(int iev=0; iev<nevt ; ++iev)
  //for(int iev=0; iev<300000 ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
    if(Centrality > cHigh || Centrality < cLow) continue;
    weight = findNcoll(Centrality) * Gen_weight;

    for(int igen = 0; igen<Gen_QQ_size; igen++){
	    JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(igen);
	    mupl_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[igen]);
	    mumi_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[igen]);

	    Double_t Rapidity_gen = fabs(JP_Gen->Rapidity());

	    //if(Rapidity_gen < 1.2) muPtCut = 3.5;
	    //else if(Rapidity_gen >1.2 && Rapidity_gen <2.1) muPtCut = (5.47-1.89*Rapidity_gen); 
	    //else if(Rapidity_gen >2.1 && Rapidity_gen <2.4) muPtCut = 1.5; 

	    if(! (JP_Gen->Pt()<50 && fabs(JP_Gen->Rapidity())<2.4 && IsAcceptanceQQ(mupl_Gen->Pt(),fabs(mupl_Gen->Eta()))&&IsAcceptanceQQ(mumi_Gen->Pt(),fabs(mumi_Gen->Eta()))) ) continue;

	    pt_weight = 1;
	    if(isPtWeight) pt_weight = fptw->Eval(JP_Gen->Pt()); 

	    if(!( fabs(JP_Gen->Rapidity())<2.4 && (mupl_Gen->Pt()>muPtCut && fabs(mupl_Gen->Eta())<2.4) && (mumi_Gen->Pt()>muPtCut && fabs(mumi_Gen->Eta())<2.4) )) continue;
	    if(Gen_mu_charge[Gen_QQ_mupl_idx[igen]]*Gen_mu_charge[Gen_QQ_mumi_idx[igen]]>0) continue;

	    if( Centrality > 0 && Centrality < 20 && JP_Gen->Pt() > 3) {
		    if(Rapidity_gen < 0.6) 					hpt_gen_1 -> Fill(JP_Gen->Pt(), weight*pt_weight);
		    else if(Rapidity_gen > 0.6 && Rapidity_gen < 1.2)	 	hpt_gen_2 -> Fill(JP_Gen->Pt(), weight*pt_weight);
		    else if(Rapidity_gen > 1.2 && Rapidity_gen < 1.6) 		hpt_gen_3 -> Fill(JP_Gen->Pt(), weight*pt_weight);
		    else if(Rapidity_gen > 1.6 && Rapidity_gen < 2.4) 		hpt_gen_4 -> Fill(JP_Gen->Pt(), weight*pt_weight);
	    }

	    if( Centrality > 20 && Centrality < 100 && JP_Gen->Pt() > 3) {
		    if(Rapidity_gen < 0.6) 					hpt_gen_5 -> Fill(JP_Gen->Pt(), weight*pt_weight);
		    else if(Rapidity_gen > 0.6 && Rapidity_gen < 1.2)	 	hpt_gen_6 -> Fill(JP_Gen->Pt(), weight*pt_weight);
		    else if(Rapidity_gen > 1.2 && Rapidity_gen < 1.6) 		hpt_gen_7 -> Fill(JP_Gen->Pt(), weight*pt_weight);
		    else if(Rapidity_gen > 1.6 && Rapidity_gen < 2.4) 		hpt_gen_8 -> Fill(JP_Gen->Pt(), weight*pt_weight);
	    }

	    if( Centrality > 100 && Centrality < 180 && JP_Gen->Pt() > 3) {
		    if(Rapidity_gen < 0.6) 					hpt_gen_9 -> Fill(JP_Gen->Pt(), weight*pt_weight);
		    else if(Rapidity_gen > 0.6 && Rapidity_gen < 1.2)	 	hpt_gen_10 -> Fill(JP_Gen->Pt(), weight*pt_weight);
		    else if(Rapidity_gen > 1.2 && Rapidity_gen < 1.6) 		hpt_gen_11 -> Fill(JP_Gen->Pt(), weight*pt_weight);
		    else if(Rapidity_gen > 1.6 && Rapidity_gen < 2.4) 		hpt_gen_12 -> Fill(JP_Gen->Pt(), weight*pt_weight);
	    }
	    if( Centrality > 0 && Centrality < 100 && JP_Gen->Pt() > 3) {
		    if(Rapidity_gen < 0.6) 					hpt_gen_13 -> Fill(JP_Gen->Pt(), weight*pt_weight);
		    else if(Rapidity_gen > 0.6 && Rapidity_gen < 1.2)	 	hpt_gen_14 -> Fill(JP_Gen->Pt(), weight*pt_weight);
		    else if(Rapidity_gen > 1.2 && Rapidity_gen < 1.6) 		hpt_gen_15 -> Fill(JP_Gen->Pt(), weight*pt_weight);
		    else if(Rapidity_gen > 1.6 && Rapidity_gen < 2.4) 		hpt_gen_16 -> Fill(JP_Gen->Pt(), weight*pt_weight);
	    }
	    if( Centrality > 0 && Centrality < 100 && JP_Gen->Pt() > 3) {
		    if(Rapidity_gen < 2.4) 					hpt_gen_17 -> Fill(JP_Gen->Pt(), weight*pt_weight);
	    }
    }


    bool HLTPass = false;
    if((HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) HLTPass=true;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
	    JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
	    mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
	    mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

      
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

      Double_t Rapidity_reco = fabs(JP_Reco->Rapidity());

      //if(Rapidity_reco < 1.2) muPtCut = 3.5;
      //else if(Rapidity_reco >1.2 && Rapidity_reco <2.1) muPtCut = (5.47-1.89*Rapidity_reco); 
      //else if(Rapidity_reco >2.1 && Rapidity_reco <2.4) muPtCut = 1.5; 

      if(! (fabs(JP_Reco->Rapidity())<2.4 && JP_Reco->Pt()<50 && IsAcceptanceQQ(mupl_Reco->Pt(),fabs(mupl_Reco->Eta()))&&IsAcceptanceQQ(mumi_Reco->Pt(),fabs(mumi_Reco->Eta()))) ) continue;

      if(!( (mupl_Reco->Pt()>muPtCut && fabs(mupl_Reco->Eta())<2.4) && (mumi_Reco->Pt()>muPtCut && fabs(mumi_Reco->Eta())<2.4) && (fabs(JP_Reco->Rapidity())<2.4 && JP_Reco->Pt()<50 && JP_Reco->Pt()>3  && JP_Reco->M()>massLow && JP_Reco->M()<massHigh))) continue;


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

       /*if( mupl_isL2 && mumi_isL3){
	 tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
	 tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);
	 SelDone = true;
	 }
	 else if( mupl_isL3 && mumi_isL2){
	 tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
	 tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
	 SelDone = true;
	 }
	 else if( mupl_isL3 && mumi_isL3){
	 int t[2] = {-1,1}; // mupl, mumi
	 int l = rand() % (2); 
       //pick up what will be L2
       if(t[l]==-1){
       tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
       tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);
       }
       else if(t[l]==1){
       tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
       tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
       }
       else {cout << "ERROR :: No random selection done !!!!" << endl; continue;}
       SelDone = true;
       }*/    
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
       tnp_weight = tnp_weight * tnp_trig_weight_mupl * tnp_trig_weight_mumi;
       if(HLTPass==true && HLTFilterPass==true){
	       counttnp++;
	       //tnp_trig_dimu = tnp_trig_weight_mupl * tnp_trig_weight_mumi;
	       tnp_trig_dimu = tnp_trig_weight;
	       hpt_tnp_trig->Fill(JP_Reco->Pt(),tnp_trig_dimu);
       }
      }

      pt_weight = 1;
      if(isPtWeight) pt_weight = fptw->Eval(JP_Reco->Pt()); 

      if(HLTPass==true && HLTFilterPass==true){
	      if( Centrality > 0 && Centrality < 20 && JP_Reco->Pt() > 3) {
	        if(Rapidity_reco < 0.6)				hpt_reco_1 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	        else if(Rapidity_reco > 0.6 && Rapidity_reco < 1.2) 	hpt_reco_2 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	        else if(Rapidity_reco > 1.2 && Rapidity_reco < 1.6) 	hpt_reco_3 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	        else if(Rapidity_reco > 1.6 && Rapidity_reco < 2.4) 	hpt_reco_4 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	      }

	      if( Centrality > 20 && Centrality < 100 && JP_Reco->Pt() > 3) {
	        if(Rapidity_reco < 0.6)				hpt_reco_5 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	        else if(Rapidity_reco > 0.6 && Rapidity_reco < 1.2) 	hpt_reco_6 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	        else if(Rapidity_reco > 1.2 && Rapidity_reco < 1.6) 	hpt_reco_7 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	        else if(Rapidity_reco > 1.6 && Rapidity_reco < 2.4) 	hpt_reco_8 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	      }

	      if( Centrality > 100 && Centrality < 180 && JP_Reco->Pt() > 3) {
	        if(Rapidity_reco < 0.6)				hpt_reco_9 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	        else if(Rapidity_reco > 0.6 && Rapidity_reco < 1.2) 	hpt_reco_10 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	        else if(Rapidity_reco > 1.2 && Rapidity_reco < 1.6) 	hpt_reco_11 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	        else if(Rapidity_reco > 1.6 && Rapidity_reco < 2.4) 	hpt_reco_12 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	      }
	      if( Centrality > 0 && Centrality < 100 && JP_Reco->Pt() > 3) {
	        if(Rapidity_reco < 0.6)					hpt_reco_13 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	        else if(Rapidity_reco > 0.6 && Rapidity_reco < 1.2) 	hpt_reco_14 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	        else if(Rapidity_reco > 1.2 && Rapidity_reco < 1.6) 	hpt_reco_15 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	        else if(Rapidity_reco > 1.6 && Rapidity_reco < 2.4) 	hpt_reco_16 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	      }
	      if( Centrality > 0 && Centrality < 100 && JP_Reco->Pt() > 3) {
	        if(Rapidity_reco < 2.4)					hpt_reco_17 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	      }
      }
      //if( Centrality > 0 && Centrality < 20 && JP_Reco->Pt() > 3) {
      //        if(Rapidity_reco < 0.6)					hpt_reco_1 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
      //        else if(Rapidity_reco > 0.6 && Rapidity_reco < 1.2) 	hpt_reco_2 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
      //        else if(Rapidity_reco > 1.2 && Rapidity_reco < 1.6) 	hpt_reco_3 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
      //        else if(Rapidity_reco > 1.6 && Rapidity_reco < 2.4) 	hpt_reco_4 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
      //}

      //if( Centrality > 20 && Centrality < 100 && JP_Reco->Pt() > 3) {
      //        if(Rapidity_reco < 0.6)					hpt_reco_5 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
      //        else if(Rapidity_reco > 0.6 && Rapidity_reco < 1.2) 	hpt_reco_6 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
      //        else if(Rapidity_reco > 1.2 && Rapidity_reco < 1.6) 	hpt_reco_7 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
      //        else if(Rapidity_reco > 1.6 && Rapidity_reco < 2.4) 	hpt_reco_8 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
      //}

      //if( Centrality > 100 && Centrality < 180 && JP_Reco->Pt() > 3) {
      //        if(Rapidity_reco < 0.6)					hpt_reco_9 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
      //        else if(Rapidity_reco > 0.6 && Rapidity_reco < 1.2) 	hpt_reco_10 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
      //        else if(Rapidity_reco > 1.2 && Rapidity_reco < 1.6) 	hpt_reco_11 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
      //        else if(Rapidity_reco > 1.6 && Rapidity_reco < 2.4) 	hpt_reco_12 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
      //}
    }
    //if(iev%100000==0) cout << ">>>>> Pt_Weight " << pt_weight << " / " << ">>>>> TnP_Weight " <<  tnp_weight << endl;
    //if(pt_weight > 0.8 && pt_weight <1.2) countPtW++;
    //if(iev%100000==0) cout << ">>>>> tnp_num " << tnp_trig_weight_num << " / " << ">>>>> tnp_den " <<  tnp_trig_weight_den <<  ">>>>> SF " << tnp_trig_weight <<endl;
    //if(iev%100000==0) cout << ">>>>> tnp_num_mL2 " << tnp_trig_weight_mumiL2_num << " / " << ">>>>> tnp_num_mL3 " << tnp_trig_weight_mumiL3_num  <<endl;
    //if(iev%100000==0) cout << ">>>>> tnp_num_pL2 " << tnp_trig_weight_muplL2_num << " / " << ">>>>> tnp_num_pL3 " << tnp_trig_weight_muplL3_num  <<endl;
    //if(iev%100000==0) cout << ">>>>> tnp_den_mL2 " << tnp_trig_weight_mumiL2_den << " / " << ">>>>> tnp_den_mL3 " << tnp_trig_weight_mumiL3_den  <<endl;
    //if(iev%100000==0) cout << ">>>>> tnp_den_pL2 " << tnp_trig_weight_muplL2_den << " / " << ">>>>> tnp_den_pL3 " << tnp_trig_weight_muplL3_den  <<endl;
    //if(iev%100000==0) cout << ">>>>> mupl_pt " << mupl_Reco->Pt() << " / " << ">>>>> mupl_eta " << mupl_Reco->Eta()  <<endl;
    //if(iev%100000==0) cout << ">>>>> mumi_pt " << mumi_Reco->Pt() << " / " << ">>>>> mumi_eta " << mumi_Reco->Eta()  <<endl;
  }

  cout << "count " << count << endl;
  cout << "counttnp " << counttnp << endl;
  //cout << "countptw " << countPtW << endl;

  //Divide
  TH1D* hpt_eff_1;
  TH1D* hpt_eff_2;
  TH1D* hpt_eff_3;
  TH1D* hpt_eff_4;

  TH1D* hpt_eff_5;
  TH1D* hpt_eff_6;
  TH1D* hpt_eff_7;
  TH1D* hpt_eff_8;

  TH1D* hpt_eff_9;
  TH1D* hpt_eff_10;
  TH1D* hpt_eff_11;
  TH1D* hpt_eff_12;

  TH1D* hpt_eff_13;
  TH1D* hpt_eff_14;
  TH1D* hpt_eff_15;
  TH1D* hpt_eff_16;
  TH1D* hpt_eff_17;

  hpt_eff_1 = (TH1D*)hpt_reco_1->Clone("hpt_eff_1");
  hpt_eff_2 = (TH1D*)hpt_reco_2->Clone("hpt_eff_2");
  hpt_eff_3 = (TH1D*)hpt_reco_3->Clone("hpt_eff_3");
  hpt_eff_4 = (TH1D*)hpt_reco_4->Clone("hpt_eff_4");

  hpt_eff_5 = (TH1D*)hpt_reco_5->Clone("hpt_eff_5");
  hpt_eff_6 = (TH1D*)hpt_reco_6->Clone("hpt_eff_6");
  hpt_eff_7 = (TH1D*)hpt_reco_7->Clone("hpt_eff_7");
  hpt_eff_8 = (TH1D*)hpt_reco_8->Clone("hpt_eff_8");
                                                 
  hpt_eff_9 = (TH1D*)hpt_reco_9->Clone("hpt_eff_9");
  hpt_eff_10 = (TH1D*)hpt_reco_10->Clone("hpt_eff_10");
  hpt_eff_11 = (TH1D*)hpt_reco_11->Clone("hpt_eff_11");
  hpt_eff_12 = (TH1D*)hpt_reco_12->Clone("hpt_eff_12");

  hpt_eff_13 = (TH1D*)hpt_reco_13->Clone("hpt_eff_13");
  hpt_eff_14 = (TH1D*)hpt_reco_14->Clone("hpt_eff_14");
  hpt_eff_15 = (TH1D*)hpt_reco_15->Clone("hpt_eff_15");
  hpt_eff_16 = (TH1D*)hpt_reco_16->Clone("hpt_eff_16");
  hpt_eff_17 = (TH1D*)hpt_reco_17->Clone("hpt_eff_17");

  hpt_eff_1->Divide(hpt_eff_1, hpt_gen_1, 1, 1, "B");
  hpt_eff_2->Divide(hpt_eff_2, hpt_gen_2, 1, 1, "B");
  hpt_eff_3->Divide(hpt_eff_3, hpt_gen_3, 1, 1, "B");
  hpt_eff_4->Divide(hpt_eff_4, hpt_gen_4, 1, 1, "B");

  hpt_eff_5->Divide(hpt_eff_5, hpt_gen_5, 1, 1, "B");
  hpt_eff_6->Divide(hpt_eff_6, hpt_gen_6, 1, 1, "B");
  hpt_eff_7->Divide(hpt_eff_7, hpt_gen_7, 1, 1, "B");
  hpt_eff_8->Divide(hpt_eff_8, hpt_gen_8, 1, 1, "B");
                                        
  hpt_eff_9->Divide(hpt_eff_9, hpt_gen_9, 1, 1, "B");
  hpt_eff_10->Divide(hpt_eff_10, hpt_gen_10, 1, 1, "B");
  hpt_eff_11->Divide(hpt_eff_11, hpt_gen_11, 1, 1, "B");
  hpt_eff_12->Divide(hpt_eff_12, hpt_gen_12, 1, 1, "B");

  hpt_eff_13->Divide(hpt_eff_13, hpt_gen_13, 1, 1, "B");
  hpt_eff_14->Divide(hpt_eff_14, hpt_gen_14, 1, 1, "B");
  hpt_eff_15->Divide(hpt_eff_15, hpt_gen_15, 1, 1, "B");
  hpt_eff_16->Divide(hpt_eff_16, hpt_gen_16, 1, 1, "B");
  hpt_eff_17->Divide(hpt_eff_17, hpt_gen_17, 1, 1, "B");


  f1->SetLineColor(kBlack);
  f1->SetLineWidth(2);
  hpt_eff_1->Fit(f1);
  hpt_eff_2->Fit(f1);
  hpt_eff_3->Fit(f1);
  hpt_eff_4->Fit(f1);

  hpt_eff_5->Fit(f1);
  hpt_eff_6->Fit(f1);
  hpt_eff_7->Fit(f1);
  hpt_eff_8->Fit(f1);
           
  hpt_eff_9->Fit(f1);
  hpt_eff_10->Fit(f1);
  hpt_eff_11->Fit(f1);
  hpt_eff_12->Fit(f1);

  hpt_eff_13->Fit(f1);
  hpt_eff_14->Fit(f1);
  hpt_eff_15->Fit(f1);
  hpt_eff_16->Fit(f1);
  hpt_eff_17->Fit(f1);

  //f1->SetLineColor(kRed+1);
  //f1->SetLineColor(kBlue+1);
  //f1->SetLineColor(kMagenta+1);
  //f1->SetLineColor(kGreen+1);
  //f1->SetLineColor(kRed-1);

  //gROOT->Macro("~/rootlogon.C");
  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.03);
  auto legend_eff_1 = new TLegend(0.6,0.84);
  auto legend_eff_2 = new TLegend(0.6,0.84);
  auto legend_eff_3 = new TLegend(0.6,0.84);
  auto legend_eff_4 = new TLegend(0.6,0.84);
  auto legend_eff_5 = new TLegend(0.6,0.84);
  auto legend_eff_6 = new TLegend(0.6,0.84);
  auto legend_eff_7 = new TLegend(0.6,0.84);
  auto legend_eff_8 = new TLegend(0.6,0.84);
  auto legend_eff_9 = new TLegend(0.6,0.84);
  auto legend_eff_10 = new TLegend(0.6,0.84);
  auto legend_eff_11 = new TLegend(0.6,0.84);
  auto legend_eff_12 = new TLegend(0.6,0.84);
  auto legend_eff_13 = new TLegend(0.6,0.84);
  auto legend_eff_14 = new TLegend(0.6,0.84);
  auto legend_eff_15 = new TLegend(0.6,0.84);
  auto legend_eff_16 = new TLegend(0.6,0.84);
  auto legend_eff_17 = new TLegend(0.6,0.84);

  gStyle->SetOptFit(0);
  //cent0-20
  TCanvas * cpt_eff_1 = new TCanvas("cpt_eff_1","cpt_eff_1",0,0,900,800);
  cpt_eff_1->cd();
  hpt_eff_1->SetTitle("Efficiency");
  hpt_eff_1->SetMarkerStyle(24);
  hpt_eff_1->SetMarkerColor(1);
  hpt_eff_1->SetLineColor(1);
  hpt_eff_1->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_1->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_1->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_1->Draw("E");
  legend_eff_1->AddEntry("hpt_eff_1","|y|: 0.0-0.6, 0-10%","lep");
  legend_eff_1->SetBorderSize(0);
  legend_eff_1->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_2 = new TCanvas("cpt_eff_2","cpt_eff_2",0,0,900,800);
  cpt_eff_2->cd();
  hpt_eff_2->SetTitle("Efficiency");
  hpt_eff_2->SetMarkerStyle(24);
  hpt_eff_2->SetMarkerColor(1);
  hpt_eff_2->SetLineColor(1);
  hpt_eff_2->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_2->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_2->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_2->Draw("E");
  legend_eff_2->AddEntry("hpt_eff_2","|y|: 0.6-1.2, 0-10%","lep");
  legend_eff_2->SetBorderSize(0);
  legend_eff_2->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_3 = new TCanvas("cpt_eff_3","cpt_eff_3",0,0,900,800);
  cpt_eff_3->cd();
  hpt_eff_3->SetTitle("Efficiency");
  hpt_eff_3->SetMarkerStyle(24);
  hpt_eff_3->SetMarkerColor(1);
  hpt_eff_3->SetLineColor(1);
  hpt_eff_3->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_3->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_3->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_3->Draw("E");
  legend_eff_3->AddEntry("hpt_eff_3","|y|: 1.2-1.6, 0-10%","lep");
  legend_eff_3->SetBorderSize(0);
  legend_eff_3->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_4 = new TCanvas("cpt_eff_4","cpt_eff_4",0,0,900,800);
  cpt_eff_4->cd();
  hpt_eff_4->SetTitle("Efficiency");
  hpt_eff_4->SetMarkerStyle(24);
  hpt_eff_4->SetMarkerColor(1);
  hpt_eff_4->SetLineColor(1);
  hpt_eff_4->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_4->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_4->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_4->Draw("E");
  legend_eff_4->AddEntry("hpt_eff_4","|y|: 1.6-2.4, 0-10%","lep");
  legend_eff_4->SetBorderSize(0);
  legend_eff_4->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  //cent20-100
  TCanvas * cpt_eff_5 = new TCanvas("cpt_eff_5","cpt_eff_5",0,0,900,800);
  cpt_eff_5->cd();
  hpt_eff_5->SetTitle("Efficiency");
  hpt_eff_5->SetMarkerStyle(24);
  hpt_eff_5->SetMarkerColor(1);
  hpt_eff_5->SetLineColor(1);
  hpt_eff_5->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_5->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_5->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_5->Draw("E");
  legend_eff_5->AddEntry("hpt_eff_5","|y|: 0.0-0.6, 10-50%","lep");
  legend_eff_5->SetBorderSize(0);
  legend_eff_5->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_6 = new TCanvas("cpt_eff_6","cpt_eff_6",0,0,900,800);
  cpt_eff_6->cd();
  hpt_eff_6->SetTitle("Efficiency");
  hpt_eff_6->SetMarkerStyle(24);
  hpt_eff_6->SetMarkerColor(1);
  hpt_eff_6->SetLineColor(1);
  hpt_eff_6->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_6->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_6->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_6->Draw("E");
  legend_eff_6->AddEntry("hpt_eff_6","|y|: 0.6-1.2, 10-50%","lep");
  legend_eff_6->SetBorderSize(0);
  legend_eff_6->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_7 = new TCanvas("cpt_eff_7","cpt_eff_7",0,0,900,800);
  cpt_eff_7->cd();
  hpt_eff_7->SetTitle("Efficiency");
  hpt_eff_7->SetMarkerStyle(24);
  hpt_eff_7->SetMarkerColor(1);
  hpt_eff_7->SetLineColor(1);
  hpt_eff_7->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_7->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_7->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_7->Draw("E");
  legend_eff_7->AddEntry("hpt_eff_7","|y|: 1.2-1.6, 10-50%","lep");
  legend_eff_7->SetBorderSize(0);
  legend_eff_7->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_8 = new TCanvas("cpt_eff_8","cpt_eff_8",0,0,900,800);
  cpt_eff_8->cd();
  hpt_eff_8->SetTitle("Efficiency");
  hpt_eff_8->SetMarkerStyle(24);
  hpt_eff_8->SetMarkerColor(1);
  hpt_eff_8->SetLineColor(1);
  hpt_eff_8->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_8->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_8->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_8->Draw("E");
  legend_eff_8->AddEntry("hpt_eff_8","|y|: 1.6-2.4, 10-50%","lep");
  legend_eff_8->SetBorderSize(0);
  legend_eff_8->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  //cent100-180
  TCanvas * cpt_eff_9 = new TCanvas("cpt_eff_9","cpt_eff_9",0,0,900,800);
  cpt_eff_9->cd();
  hpt_eff_9->SetTitle("Efficiency");
  hpt_eff_9->SetMarkerStyle(24);
  hpt_eff_9->SetMarkerColor(1);
  hpt_eff_9->SetLineColor(1);
  hpt_eff_9->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_9->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_9->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_9->Draw("E");
  legend_eff_9->AddEntry("hpt_eff_9","|y|: 0.0-0.6, 50-90%","lep");
  legend_eff_9->SetBorderSize(0);
  legend_eff_9->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_10 = new TCanvas("cpt_eff_10","cpt_eff_10",0,0,900,800);
  cpt_eff_10->cd();
  hpt_eff_10->SetTitle("Efficiency");
  hpt_eff_10->SetMarkerStyle(24);
  hpt_eff_10->SetMarkerColor(1);
  hpt_eff_10->SetLineColor(1);
  hpt_eff_10->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_10->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_10->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_10->Draw("E");
  legend_eff_10->AddEntry("hpt_eff_10","|y|: 0.6-1.2, 50-90%","lep");
  legend_eff_10->SetBorderSize(0);
  legend_eff_10->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_11 = new TCanvas("cpt_eff_11","cpt_eff_11",0,0,900,800);
  cpt_eff_11->cd();
  hpt_eff_11->SetTitle("Efficiency");
  hpt_eff_11->SetMarkerStyle(24);
  hpt_eff_11->SetMarkerColor(1);
  hpt_eff_11->SetLineColor(1);
  hpt_eff_11->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_11->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_11->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_11->Draw("E");
  legend_eff_11->AddEntry("hpt_eff_11","|y|: 1.2-1.6, 50-90%","lep");
  legend_eff_11->SetBorderSize(0);
  legend_eff_11->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_12 = new TCanvas("cpt_eff_12","cpt_eff_12",0,0,900,800);
  cpt_eff_12->cd();
  hpt_eff_12->SetTitle("Efficiency");
  hpt_eff_12->SetMarkerStyle(24);
  hpt_eff_12->SetMarkerColor(1);
  hpt_eff_12->SetLineColor(1);
  hpt_eff_12->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_12->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_12->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_12->Draw("E");
  legend_eff_12->AddEntry("hpt_eff_12","|y|: 1.6-2.4, 50-90%","lep");
  legend_eff_12->SetBorderSize(0);
  legend_eff_12->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  //cent0-100
  TCanvas * cpt_eff_13 = new TCanvas("cpt_eff_13","cpt_eff_13",0,0,900,800);
  cpt_eff_13->cd();
  hpt_eff_13->SetTitle("Efficiency");
  hpt_eff_13->SetMarkerStyle(24);
  hpt_eff_13->SetMarkerColor(1);
  hpt_eff_13->SetLineColor(1);
  hpt_eff_13->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_13->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_13->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_13->Draw("E");
  legend_eff_13->AddEntry("hpt_eff_13","|y|: 0.0-0.6, 0-50%","lep");
  legend_eff_13->SetBorderSize(0);
  legend_eff_13->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_14 = new TCanvas("cpt_eff_14","cpt_eff_14",0,0,900,800);
  cpt_eff_14->cd();
  hpt_eff_14->SetTitle("Efficiency");
  hpt_eff_14->SetMarkerStyle(24);
  hpt_eff_14->SetMarkerColor(1);
  hpt_eff_14->SetLineColor(1);
  hpt_eff_14->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_14->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_14->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_14->Draw("E");
  legend_eff_14->AddEntry("hpt_eff_14","|y|: 0.6-1.2, 0-50%","lep");
  legend_eff_14->SetBorderSize(0);
  legend_eff_14->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_15 = new TCanvas("cpt_eff_15","cpt_eff_15",0,0,900,800);
  cpt_eff_15->cd();
  hpt_eff_15->SetTitle("Efficiency");
  hpt_eff_15->SetMarkerStyle(24);
  hpt_eff_15->SetMarkerColor(1);
  hpt_eff_15->SetLineColor(1);
  hpt_eff_15->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_15->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_15->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_15->Draw("E");
  legend_eff_15->AddEntry("hpt_eff_15","|y|: 1.2-1.6, 0-50%","lep");
  legend_eff_15->SetBorderSize(0);
  legend_eff_15->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_16 = new TCanvas("cpt_eff_16","cpt_eff_16",0,0,900,800);
  cpt_eff_16->cd();
  hpt_eff_16->SetTitle("Efficiency");
  hpt_eff_16->SetMarkerStyle(24);
  hpt_eff_16->SetMarkerColor(1);
  hpt_eff_16->SetLineColor(1);
  hpt_eff_16->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_16->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_16->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_16->Draw("E");
  legend_eff_16->AddEntry("hpt_eff_16","|y|: 1.6-2.4, 0-50%","lep");
  legend_eff_16->SetBorderSize(0);
  legend_eff_16->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");

  TCanvas * cpt_eff_17 = new TCanvas("cpt_eff_17","cpt_eff_17",0,0,900,800);
  cpt_eff_17->cd();
  hpt_eff_17->SetTitle("Efficiency");
  hpt_eff_17->SetMarkerStyle(24);
  hpt_eff_17->SetMarkerColor(1);
  hpt_eff_17->SetLineColor(1);
  hpt_eff_17->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hpt_eff_17->GetYaxis()->SetTitle("Efficiency");
  hpt_eff_17->GetYaxis()->SetRangeUser(0.,1.3);
  hpt_eff_17->Draw("E");
  legend_eff_17->AddEntry("hpt_eff_17","|y|: 0-2.4, 0-50%","lep");
  legend_eff_17->SetBorderSize(0);
  legend_eff_17->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");


  hpt_eff_1 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to20_absy0_0p6  ",isTnP, isPtWeight));
  hpt_eff_2 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to20_absy0p6_1p2",isTnP, isPtWeight));
  hpt_eff_3 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to20_absy1p2_1p6",isTnP, isPtWeight));
  hpt_eff_4 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to20_absy1p6_2p4",isTnP, isPtWeight));
  hpt_eff_5 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent20to100_absy0_0p6  ",isTnP, isPtWeight));
  hpt_eff_6 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent20to100_absy0p6_1p2",isTnP, isPtWeight));
  hpt_eff_7 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent20to100_absy1p2_1p6",isTnP, isPtWeight));
  hpt_eff_8 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent20to100_absy1p6_2p4",isTnP, isPtWeight));
  hpt_eff_9 ->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent100to180_absy0_0p6  ",isTnP, isPtWeight));
  hpt_eff_10->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent100to180_absy0p6_1p2",isTnP, isPtWeight));
  hpt_eff_11->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent100to180_absy1p2_1p6",isTnP, isPtWeight));
  hpt_eff_12->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent100to180_absy1p6_2p4",isTnP, isPtWeight));
  hpt_eff_13->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to100_absy0_0p6  ",isTnP, isPtWeight));
  hpt_eff_14->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to100_absy0p6_1p2",isTnP, isPtWeight));
  hpt_eff_15->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to100_absy1p2_1p6",isTnP, isPtWeight));
  hpt_eff_16->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to100_absy1p6_2p4",isTnP, isPtWeight));
  hpt_eff_17->SetName(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent0to100_absy0_2p4",isTnP, isPtWeight));

  TString outFileName = Form("mc_eff_vs_pt_cent_rap_prompt_pbpb_Jpsi_PtW%d_tnp%d.root",isPtWeight,isTnP); 
  if(state==2) outFileName = Form("mc_eff_vs_pt_cent_rap_nprompt_pbpb_Jpsi_PtW%d_tnp%d.root",isPtWeight,isTnP); 
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hpt_eff_1->Write();
  hpt_eff_2->Write();
  hpt_eff_3->Write();
  hpt_eff_4->Write();
  hpt_eff_5->Write();
  hpt_eff_6->Write();
  hpt_eff_7->Write();
  hpt_eff_8->Write();
  hpt_eff_9->Write();
  hpt_eff_10->Write();
  hpt_eff_11->Write();
  hpt_eff_12->Write();
  hpt_eff_13->Write();
  hpt_eff_14->Write();
  hpt_eff_15->Write();
  hpt_eff_16->Write();
  hpt_eff_17->Write();
  f1->Write();

  if(isTnP) hpt_tnp_trig->Write();
  outFile->Close();

}
