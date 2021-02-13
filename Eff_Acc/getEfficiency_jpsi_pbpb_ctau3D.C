#include <iostream>

#include <TLorentzVector.h>
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../HiEvtPlaneList.h"
#include "../cutsAndBin.h"
#include "../Style_jaebeom.h"
#include "tnp_weight_lowptPbPb.h"
#include <TAttMarker.h>

using namespace std;

void getEfficiency_jpsi_pbpb_ctau3D(
  float ptLow = 0.0, float ptHigh = 50.0,
  float yLow = 1.6, float yHigh = 2.4,
  int cLow = 0, int cHigh = 180, bool isTnP = true, bool isPtWeight = false, int state=1
  ) {

  gStyle->SetOptStat(0);
  int kTrigSel = 12;	//jpsi=12,Upsilon=13

  float muPtCut = 0; //3.5, 1.8
  float muEtaCut = 2.4;

  ////Upsilon
  //float massLow = 8.0;
  //float massHigh = 10.0;
  //if(state==2){massLow = 8.5; massHigh = 11;}

  //jpsi
  float massLow = 2.6;
  float massHigh = 3.5;
  if(state==2){massLow = 2.6; massHigh = 3.5;}


  double min = 0;
  double max = ptHigh;
  double binwidth = 1;
  //const int numBins = 31; //50;//(max-min)/binwidth;  //31//
  const int numBins = 6; //50;//(max-min)/binwidth;  //31//

  //input files
  //PbPb
  TString inputMC = "/work2/Oniatree/JPsi/OniatreeMC_JPsi_pThat2_TunedCP5_HydjetDrumMB_5p02TeV_Pythia8_part*.root";	//PbPb_prompt
  if(state==2) inputMC = "/work2/Oniatree/JPsi/OniatreeMC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8.root";	//PbPb_non prompt
  TChain* mytree = new TChain("myTree"); 
  mytree->Add(inputMC.Data());
  TString outFileName = "mc_eff_vs_ctau3D_cent_prompt_pbpb_Jpsi.root"; 
  if(state==2) outFileName = "mc_eff_vs_ctau3D_cent_nprompt_pbpb_Jpsi.root"; 


  //pT reweighting function
  //TFile *fPtW = new TFile(Form("../Reweight/WeightedFunc/Func_dNdpT_%dS.root",state),"read");
  //TF1* f1 = (TF1*) fPtW->Get("fitRatio");

  
  //double ptBin[numBins+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27,30,34,38,42,46,50};
  double ptBin[numBins+1] = {0,3,4.5,6.5,9,12,50};
  float yBin[7] = {0,0.4,0.8,1.2,1.6,2.0,2.4};

  TH1D* hctau3D_reco_1 = new TH1D("hctau3D_reco_1","hctau3D_reco_1",10,0,10);
  TH1D* hctau3D_reco_2 = new TH1D("hctau3D_reco_2","hctau3D_reco_2",10,0,10);
  TH1D* hctau3D_reco_3 = new TH1D("hctau3D_reco_3","hctau3D_reco_3",10,0,10);
  TH1D* hctau3D_reco_4 = new TH1D("hctau3D_reco_4","hctau3D_reco_4",10,0,10);
  TH1D* hctau3D_reco_5 = new TH1D("hctau3D_reco_5","hctau3D_reco_5",10,0,10);

  TH1D* hctau3D_gen_1 = new TH1D("hctau3D_gen_1","hctau3D_gen_1",10,0,10);
  TH1D* hctau3D_gen_2 = new TH1D("hctau3D_gen_2","hctau3D_gen_2",10,0,10);
  TH1D* hctau3D_gen_3 = new TH1D("hctau3D_gen_3","hctau3D_gen_3",10,0,10);
  TH1D* hctau3D_gen_4 = new TH1D("hctau3D_gen_4","hctau3D_gen_4",10,0,10);
  TH1D* hctau3D_gen_5 = new TH1D("hctau3D_gen_5","hctau3D_gen_5",10,0,10);

  TH1D* hy_gen = new TH1D("hy_gen","hy_gen",6,yBin);
  TH1D* hy_reco = new TH1D("hy_reco","hy_reco",6,yBin);
  hy_gen ->Sumw2();
  hy_reco ->Sumw2();


  hctau3D_reco_1->Sumw2();
  hctau3D_reco_2->Sumw2();
  hctau3D_reco_3->Sumw2();
  hctau3D_reco_4->Sumw2();
  hctau3D_reco_5->Sumw2();

  hctau3D_gen_1->Sumw2();
  hctau3D_gen_2->Sumw2();
  hctau3D_gen_3->Sumw2();
  hctau3D_gen_4->Sumw2();
  hctau3D_gen_5->Sumw2();

  //hctau3D_reco_1->SetTitle("Reco: Rapidity 0.0-2.4");
  //hctau3D_reco_2->SetTitle("Reco: Rapidity 0.0-0.6");
  //hctau3D_reco_3->SetTitle("Reco: Rapidity 0.6-1.2");
  //hctau3D_reco_4->SetTitle("Reco: Rapidity 1.2-1.8");
  //hctau3D_reco_5->SetTitle("Reco: Rapidity 1.8-2.4");

  hctau3D_reco_1->SetTitle("Reco: Rapidity 1.6-2.4 Pt 3-4.5");
  hctau3D_reco_2->SetTitle("Reco: Rapidity 1.6-2.4 Pt 4.5-6.5");
  hctau3D_reco_3->SetTitle("Reco: Rapidity 1.6-2.4 Pt 6.5-9");
  hctau3D_reco_4->SetTitle("Reco: Rapidity 1.6-2.4 Pt 9-12");
  hctau3D_reco_5->SetTitle("Reco: Rapidity 1.6-2.4 Pt 12-50");

  hctau3D_reco_1->GetXaxis()->SetTitle("ctau3D");
  hctau3D_reco_2->GetXaxis()->SetTitle("ctau3D");
  hctau3D_reco_3->GetXaxis()->SetTitle("ctau3D");
  hctau3D_reco_4->GetXaxis()->SetTitle("ctau3D");
  hctau3D_reco_5->GetXaxis()->SetTitle("ctau3D");

  hctau3D_gen_1->SetTitle("Gen: Rapidity 1.6-2.4 Pt 3-4.5");
  hctau3D_gen_2->SetTitle("Gen: Rapidity 1.6-2.4 Pt 4.5-6.5");
  hctau3D_gen_3->SetTitle("Gen: Rapidity 1.6-2.4 Pt 6.5-9");
  hctau3D_gen_4->SetTitle("Gen: Rapidity 1.6-2.4 Pt 9-12");
  hctau3D_gen_5->SetTitle("Gen: Rapidity 1.6-2.4 Pt 12-50");

  hctau3D_gen_1->GetXaxis()->SetTitle("ctau3D");
  hctau3D_gen_2->GetXaxis()->SetTitle("ctau3D");
  hctau3D_gen_3->GetXaxis()->SetTitle("ctau3D");
  hctau3D_gen_4->GetXaxis()->SetTitle("ctau3D");
  hctau3D_gen_5->GetXaxis()->SetTitle("ctau3D");

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
  Float_t		 Gen_QQ_ctau3D;   //!
  TBranch        *b_Gen_QQ_ctau3D;   //!

  Gen_QQ_4mom = 0; Gen_mu_4mom = 0;
  mytree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  mytree->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
  mytree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  mytree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  mytree->SetBranchAddress("Gen_QQ_ctau3D", &Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);

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
  Float_t		 Reco_QQ_ctau3D;   //!
  TBranch        *b_Reco_QQ_ctau3D;   //!

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
  mytree->SetBranchAddress("Reco_QQ_ctau3D", &Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);

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
  double tnp_trig_weight_mupl = -1;
  double tnp_trig_weight_mumi = -1;
  double pt_weight = 1;
  
  double tnp_trig_dimu=-1;
  TH2D* hpt_tnp_trig = new TH2D("hpt_tnp_trig","hpt_tnp_trig",numBins,ptBin,40,0,2);

  int kL2filter = 16;	//jpsi=16,Upsilon=38
  int kL3filter = 17;	//jpsi=17,Upsilon=39

  int count =0;
  int counttnp =0;
  const int nevt = mytree->GetEntries();
  cout << "Total Events : " << nevt << endl;
  for(int iev=0; iev<nevt ; ++iev)
  //cout << "Total Events : " << 500000 << endl;
  //for(int iev=0; iev<500000 ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;
    //if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << 500000 <<  " ("<<(int)(100.*iev/500000) << "%)" << endl;

    mytree->GetEntry(iev);
    if(Centrality > cHigh || Centrality < cLow) continue;
    weight = findNcoll(Centrality) * Gen_weight;

    for(int igen = 0; igen<Gen_QQ_size; igen++){
	    JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(igen);
	    mupl_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[igen]);
	    mumi_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[igen]);

	    Double_t Rapidity_g = fabs(JP_Gen->Rapidity());

	    if(Rapidity_g < 1.2) muPtCut = 3.5;
	    else if(Rapidity_g >1.2 && Rapidity_g <2.1) muPtCut = (5.77-1.89*Rapidity_g); 
	    else if(Rapidity_g >2.1 && Rapidity_g <2.4) muPtCut = 1.8; 

	    if(!( fabs(JP_Gen->Rapidity())<2.4 && (mupl_Gen->Pt()>muPtCut && fabs(mupl_Gen->Eta())<2.4) && (mumi_Gen->Pt()>muPtCut && fabs(mumi_Gen->Eta())<2.4) )) continue;
	    if(Gen_mu_charge[Gen_QQ_mupl_idx[igen]]*Gen_mu_charge[Gen_QQ_mumi_idx[igen]]>0) continue;

	    pt_weight = 1;
	    //if(isPtWeight) pt_weight = f1->Eval(JP_Gen->Pt()); 

	    hy_gen->Fill(Rapidity_g, weight);
	    //hctau3D_gen_1->Fill(JP_Gen->Pt(),weight*pt_weight);
	    //if(Centrality < 20) hctau3D_gen_2 -> Fill(JP_Gen->Pt(), weight*pt_weight);
	    //else if(Centrality > 20 && Centrality < 60) hctau3D_gen_3 -> Fill(JP_Gen->Pt(), weight*pt_weight);
	    //else if(Centrality > 60 && Centrality < 120) hctau3D_gen_4 -> Fill(JP_Gen->Pt(), weight*pt_weight);
	    //else if(Centrality > 120 && Centrality < 180) hctau3D_gen_5 -> Fill(JP_Gen->Pt(), weight*pt_weight);
	    //
	    if(JP_Gen->Pt()>3 && JP_Gen->Pt()<4.5)hctau3D_gen_1->Fill(Gen_QQ_ctau3D,weight*pt_weight);
	    else if(JP_Gen->Pt()>4.5 && JP_Gen->Pt()<6.5) hctau3D_gen_2 -> Fill(Gen_QQ_ctau3D, weight*pt_weight);
	    else if(JP_Gen->Pt()>6.5 && JP_Gen->Pt()<9) hctau3D_gen_3 -> Fill(Gen_QQ_ctau3D, weight*pt_weight);
	    else if(JP_Gen->Pt()>9 && JP_Gen->Pt()<12) hctau3D_gen_4 -> Fill(Gen_QQ_ctau3D, weight*pt_weight);
	    else if(JP_Gen->Pt()>12 && JP_Gen->Pt()<50) hctau3D_gen_5 -> Fill(Gen_QQ_ctau3D, weight*pt_weight);
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

      Double_t Rapidity = fabs(JP_Reco->Rapidity());

      if(Rapidity < 1.2) muPtCut = 3.5;
      else if(Rapidity >1.2 && Rapidity <2.1) muPtCut = (5.77-1.89*Rapidity); 
      else if(Rapidity >2.1 && Rapidity <2.4) muPtCut = 1.8; 

      if(!( (mupl_Reco->Pt()>muPtCut && fabs(mupl_Reco->Eta())<2.4) && (mumi_Reco->Pt()>muPtCut && fabs(mumi_Reco->Eta())<2.4) && (fabs(JP_Reco->Rapidity())<2.4 && JP_Reco->Pt()<50  && JP_Reco->M()>massLow && JP_Reco->M()<massHigh))) continue;

      //hy_reco->Fill(Rapidity, weight*pt_weight);
      //if(JP_Reco->Pt()>3 && JP_Reco->Pt()<4.5)hctau3D_reco_1->Fill(Reco_QQ_ctau3D,weight*pt_weight);
      //else if(JP_Reco->Pt()>4.5 && JP_Reco->Pt()<6.5) hctau3D_reco_2->Fill(Reco_QQ_ctau3D, weight*pt_weight);
      //else if(JP_Reco->Pt()>6.5 && JP_Reco->Pt()<9) hctau3D_reco_3->Fill(Reco_QQ_ctau3D, weight*pt_weight);
      //else if(JP_Reco->Pt()>9 && JP_Reco->Pt()<12) hctau3D_reco_4->Fill(Reco_QQ_ctau3D, weight*pt_weight);
      //else if(JP_Reco->Pt()>12 && JP_Reco->Pt()<50) hctau3D_reco_5->Fill(Reco_QQ_ctau3D, weight*pt_weight);
      //hctau3D_reco_1->Fill(JP_Reco->Pt(),weight * pt_weight);
      //if(Rapidity < 0.6) hctau3D_reco_2 -> Fill(JP_Reco->Pt(), weight *pt_weight);
      //else if(Rapidity > 0.6 && Rapidity < 1.2) hctau3D_reco_3 -> Fill(JP_Reco->Pt(), weight *pt_weight);
      //else if(Rapidity > 1.2 && Rapidity < 1.8) hctau3D_reco_4 -> Fill(JP_Reco->Pt(), weight *pt_weight);
      //else if(Rapidity > 1.8 && Rapidity < 2.4) hctau3D_reco_5 -> Fill(JP_Reco->Pt(), weight *pt_weight);
      if(HLTPass==true && HLTFilterPass==true) count++;
      if(isTnP){
       tnp_weight = 1;
       tnp_trig_weight_mupl = -1;
       tnp_trig_weight_mumi = -1;
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

       if( mupl_isL2 && mumi_isL3){
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
       }

       //if(SelDone == false){cout << "ERROR :: No muon filter combination selected !!!!" << endl; continue;}
       //if((tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1)){cout << "ERROR :: No trigger muon tnp scale factors selected !!!!" << endl; continue;}
       tnp_weight = tnp_weight * tnp_trig_weight_mupl * tnp_trig_weight_mumi;
       if(HLTPass==true && HLTFilterPass==true){
         counttnp++;
	 tnp_trig_dimu = tnp_trig_weight_mupl * tnp_trig_weight_mumi;
	 hpt_tnp_trig->Fill(JP_Reco->Pt(),tnp_trig_dimu);
       }
      }

      pt_weight = 1;
      //if(isPtWeight) pt_weight = f1->Eval(JP_Reco->Pt());

      if(HLTPass==true && HLTFilterPass==true){
	      //hy_reco->Fill(Rapidity, weight);
	      //hpt_reco_Trig_1->Fill(JP_Reco->Pt(),weight* tnp_weight* pt_weight);
	      //if(Centrality < 20) hpt_reco_Trig_2 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	      //else if(Centrality > 20 && Centrality < 60) hpt_reco_Trig_3 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	      //else if(Centrality > 60 && Centrality < 120) hpt_reco_Trig_4 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	      ////else if(Centrality > 120 && Centrality < 180) hpt_reco_Trig_5 -> Fill(JP_Reco->Pt(), weight* tnp_weight *pt_weight);
	      //else if(Rapidity > 1.6 && Rapidity < 2.4) hpt_reco_Trig_5 -> Fill(JP_Reco->Pt(), weight * tnp_weight *pt_weight);
	      hy_reco->Fill(Rapidity, weight* tnp_weight* pt_weight);
	      if(JP_Reco->Pt()>3 && JP_Reco->Pt()<4.5)hctau3D_reco_1->Fill(Reco_QQ_ctau3D,weight*tnp_weight*pt_weight);
	      else if(JP_Reco->Pt()>4.5 && JP_Reco->Pt()<6.5) hctau3D_reco_2->Fill(Reco_QQ_ctau3D, weight*tnp_weight*pt_weight);
	      else if(JP_Reco->Pt()>6.5 && JP_Reco->Pt()<9) hctau3D_reco_3->Fill(Reco_QQ_ctau3D, weight*tnp_weight*pt_weight);
	      else if(JP_Reco->Pt()>9 && JP_Reco->Pt()<12) hctau3D_reco_4->Fill(Reco_QQ_ctau3D, weight*tnp_weight*pt_weight);
	      else if(JP_Reco->Pt()>12 && JP_Reco->Pt()<50) hctau3D_reco_5->Fill(Reco_QQ_ctau3D, weight*tnp_weight*pt_weight);
      }

    }
  }

  cout << "count " << count << endl;
  cout << "counttnp " << counttnp << endl;
  


  //Draw
  TCanvas * cctau3D_reco = new TCanvas("cctau3D_reco","cctau3D_reco",0,0,1000,1000);
  cctau3D_reco->Divide(2,3);
  cctau3D_reco->cd(1);
  hctau3D_reco_1->Draw();
  cctau3D_reco->cd(2);
  hctau3D_reco_2->Draw();
  cctau3D_reco->cd(3);
  hctau3D_reco_3->Draw();
  cctau3D_reco->cd(4);
  hctau3D_reco_4->Draw();
  cctau3D_reco->cd(5);
  hctau3D_reco_5->Draw();

  //TCanvas * cctau3D_reco_1 = new TCanvas("cctau3D_reco_1","cctau3D_reco_1",0,0,400,400);
  //cctau3D_reco_1->cd();
  //hctau3D_reco_1->Draw();

  //TCanvas * cctau3D_reco_2 = new TCanvas("cctau3D_reco_2","cctau3D_reco_2",0,0,400,400);
  //cctau3D_reco_2->cd();
  //hctau3D_reco_2->Draw();

  //TCanvas * cctau3D_reco_3 = new TCanvas("cctau3D_reco_3","cctau3D_reco_3",0,0,400,400);
  //cctau3D_reco_3->cd();
  //hctau3D_reco_3->Draw();

  //TCanvas * cctau3D_reco_4 = new TCanvas("cctau3D_reco_4","cctau3D_reco_4",0,0,400,400);
  //cctau3D_reco_4->cd();
  //hctau3D_reco_4->Draw();

  //TCanvas * cctau3D_reco_5 = new TCanvas("cctau3D_reco_5","cctau3D_reco_5",0,0,400,400);
  //cctau3D_reco_5->cd();
  //hctau3D_reco_5->Draw();

  //Gen
  //TCanvas * cy_gen = new TCanvas("cctau3D_gen","cctau3D_gen",0,0,400,400);
  //cy_gen->cd();
  //hy_gen->Draw();

  TCanvas * cctau3D_gen = new TCanvas("cctau3D_gen","cctau3D_gen",0,0,1000,1000);
  cctau3D_gen->Divide(2,3);
  cctau3D_gen->cd(1);
  hctau3D_gen_1->Draw();
  cctau3D_gen->cd(2);
  hctau3D_gen_2->Draw();
  cctau3D_gen->cd(3);
  hctau3D_gen_3->Draw();
  cctau3D_gen->cd(4);
  hctau3D_gen_4->Draw();
  cctau3D_gen->cd(5);
  hctau3D_gen_5->Draw();

  //TCanvas * cctau3D_gen_1 = new TCanvas("cctau3D_gen_1","cctau3D_gen_1",0,0,400,400);
  //cctau3D_gen_1->cd();
  //hctau3D_gen_1->Draw();

  //TCanvas * cctau3D_gen_2 = new TCanvas("cctau3D_gen_2","cctau3D_gen_2",0,0,400,400);
  //cctau3D_gen_2->cd();
  //hctau3D_gen_2->Draw();

  //TCanvas * cctau3D_gen_3 = new TCanvas("cctau3D_gen_3","cctau3D_gen_3",0,0,400,400);
  //cctau3D_gen_3->cd();
  //hctau3D_gen_3->Draw();

  //TCanvas * cctau3D_gen_4 = new TCanvas("cctau3D_gen_4","cctau3D_gen_4",0,0,400,400);
  //cctau3D_gen_4->cd();
  //hctau3D_gen_4->Draw();

  //TCanvas * cctau3D_gen_5 = new TCanvas("cctau3D_gen_5","cctau3D_gen_5",0,0,400,400);
  //cctau3D_gen_5->cd();
  //hctau3D_gen_5->Draw();

  //cctau3D_reco_1->Update();
  //cctau3D_reco_2->Update();
  //cctau3D_reco_3->Update();
  //cctau3D_reco_4->Update();
  //cctau3D_reco_5->Update();

  //Divide
  TH1D* hctau3D_eff_1;
  TH1D* hctau3D_eff_2;
  TH1D* hctau3D_eff_3;
  TH1D* hctau3D_eff_4;
  TH1D* hctau3D_eff_5;

  hctau3D_eff_1 = (TH1D*)hctau3D_reco_1->Clone("hctau3D_eff_1");
  hctau3D_eff_2 = (TH1D*)hctau3D_reco_2->Clone("hctau3D_eff_2");
  hctau3D_eff_3 = (TH1D*)hctau3D_reco_3->Clone("hctau3D_eff_3");
  hctau3D_eff_4 = (TH1D*)hctau3D_reco_4->Clone("hctau3D_eff_4");
  hctau3D_eff_5 = (TH1D*)hctau3D_reco_5->Clone("hctau3D_eff_5");

  hctau3D_eff_1->Divide(hctau3D_eff_1, hctau3D_gen_1, 1, 1, "B");
  hctau3D_eff_2->Divide(hctau3D_eff_2, hctau3D_gen_2, 1, 1, "B");
  hctau3D_eff_3->Divide(hctau3D_eff_3, hctau3D_gen_3, 1, 1, "B");
  hctau3D_eff_4->Divide(hctau3D_eff_4, hctau3D_gen_4, 1, 1, "B");
  hctau3D_eff_5->Divide(hctau3D_eff_5, hctau3D_gen_5, 1, 1, "B");

  TH1D* hy_eff;
  hy_eff = (TH1D*)hy_reco->Clone("hy_eff");
  hy_eff->Divide(hy_eff, hy_gen, 1, 1, "B");

  hctau3D_eff_1->SetTitle("Eff: Rapidity 1.6-2.4 Pt 3-4.5");
  hctau3D_eff_2->SetTitle("Eff: Rapidity 1.6-2.4 Pt 4.5-6.5");
  hctau3D_eff_3->SetTitle("Eff: Rapidity 1.6-2.4 Pt 6.5-9");
  hctau3D_eff_4->SetTitle("Eff: Rapidity 1.6-2.4 Pt 9-12");
  hctau3D_eff_5->SetTitle("Eff: Rapidity 1.6-2.4 Pt 12-50");

  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.03);
  auto legend1 = new TLegend(0.6,0.84);
  auto legend2 = new TLegend(0.6,0.84);
  auto legend3 = new TLegend(0.6,0.84);
  auto legend4 = new TLegend(0.6,0.84);
  auto legend5 = new TLegend(0.6,0.84);

  TCanvas * cctau3D_eff = new TCanvas("cctau3D_eff","cctau3D_eff",0,0,1000,1000);
  cctau3D_eff->Divide(2,3);
  cctau3D_eff->cd(1);
  hctau3D_eff_1->SetMarkerStyle(24);
  hctau3D_eff_1->SetMarkerColor(1);
  hctau3D_eff_1->SetLineColor(1);
  hctau3D_eff_1->GetXaxis()->SetTitle("ctau3D");
  hctau3D_eff_1->GetYaxis()->SetTitle("Efficiency");
  hctau3D_eff_1->GetYaxis()->SetRangeUser(0.,1.1);
  hctau3D_eff_1->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");
  legend1->AddEntry("hctau3D_eff_1","|y|:1.6-2.4 Pt:3-4.5","lep");
  legend1->SetBorderSize(0);
  legend1->Draw("same");

  cctau3D_eff->cd(2);
  hctau3D_eff_2->SetMarkerStyle(24);
  hctau3D_eff_2->SetMarkerColor(1);
  hctau3D_eff_2->SetLineColor(1);
  hctau3D_eff_2->GetXaxis()->SetTitle("ctau3D");
  hctau3D_eff_2->GetYaxis()->SetTitle("Efficiency");
  hctau3D_eff_2->GetYaxis()->SetRangeUser(0.,1.1);
  hctau3D_eff_2->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");
  legend2->AddEntry("hctau3D_eff_1","|y|:1.6-2.4 Pt:4.5-6.5","lep");
  legend2->SetBorderSize(0);
  legend2->Draw("same");

  cctau3D_eff->cd(3);
  hctau3D_eff_3->SetMarkerStyle(24);
  hctau3D_eff_3->SetMarkerColor(1);
  hctau3D_eff_3->SetLineColor(1);
  hctau3D_eff_3->GetXaxis()->SetTitle("ctau3D");
  hctau3D_eff_3->GetYaxis()->SetTitle("Efficiency");
  hctau3D_eff_3->GetYaxis()->SetRangeUser(0.,1.1);
  hctau3D_eff_3->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");
  legend3->AddEntry("hctau3D_eff_1","|y|:1.6-2.4 Pt:6.5-9","lep");
  legend3->SetBorderSize(0);
  legend3->Draw("same");

  cctau3D_eff->cd(4);
  hctau3D_eff_4->SetMarkerStyle(24);
  hctau3D_eff_4->SetMarkerColor(1);
  hctau3D_eff_4->SetLineColor(1);
  hctau3D_eff_4->GetXaxis()->SetTitle("ctau3D");
  hctau3D_eff_4->GetYaxis()->SetTitle("Efficiency");
  hctau3D_eff_4->GetYaxis()->SetRangeUser(0.,1.1);
  hctau3D_eff_4->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");
  legend4->AddEntry("hctau3D_eff_1","|y|:1.6-2.4 Pt:9-12","lep");
  legend4->SetBorderSize(0);
  legend4->Draw("same");

  cctau3D_eff->cd(5);
  hctau3D_eff_5->SetMarkerStyle(24);
  hctau3D_eff_5->SetMarkerColor(1);
  hctau3D_eff_5->SetLineColor(1);
  hctau3D_eff_5->GetXaxis()->SetTitle("ctau3D");
  hctau3D_eff_5->GetYaxis()->SetTitle("Efficiency");
  hctau3D_eff_5->GetYaxis()->SetRangeUser(0.,1.1);
  hctau3D_eff_5->Draw("E");
  lt1->SetTextSize(0.03);
  if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");
  legend5->AddEntry("hctau3D_eff_1","|y|:1.6-2.4 Pt:12-50","lep");
  legend5->SetBorderSize(0);
  legend5->Draw("same");


  //TCanvas * cctau3D_eff_1 = new TCanvas("cctau3D_eff_1","cctau3D_eff_1",0,0,400,400);
  //cctau3D_eff_1->cd();
  //hctau3D_eff_1->Draw();

  //TCanvas * cctau3D_eff_2 = new TCanvas("cctau3D_eff_2","cctau3D_eff_2",0,0,400,400);
  //cctau3D_eff_2->cd();
  //hctau3D_eff_2->Draw();

  //TCanvas * cctau3D_eff_3 = new TCanvas("cctau3D_eff_3","cctau3D_eff_3",0,0,400,400);
  //cctau3D_eff_3->cd();
  //hctau3D_eff_3->Draw();

  //TCanvas * cctau3D_eff_4 = new TCanvas("cctau3D_eff_4","cctau3D_eff_4",0,0,400,400);
  //cctau3D_eff_4->cd();
  //hctau3D_eff_4->Draw();

  //TCanvas * cctau3D_eff_5 = new TCanvas("cctau3D_eff_5","cctau3D_eff_5",0,0,400,400);
  //cctau3D_eff_5->cd();
  //hctau3D_eff_5->Draw();
//draw same
  ////gROOT->Macro("~/rootlogon.C");
  //TLatex *lt1 = new TLatex();
  //lt1->SetNDC();
  //lt1->SetTextSize(0.03);
  //auto legend = new TLegend(0.6,0.84);

  //gStyle->SetOptFit(0);
  //TCanvas * cctau3D_eff = new TCanvas("cctau3D_eff","cctau3D_eff",0,0,900,800);
  //cctau3D_eff->cd();
  //hctau3D_eff_1->SetMarkerStyle(24);
  ////hctau3D_eff_1->SetMarkerSize(7);
  //hctau3D_eff_1->SetMarkerColor(1);
  //hctau3D_eff_1->SetLineColor(1);
  ////hctau3D_eff_1->GetXaxis()->CenterTitle();
  //hctau3D_eff_1->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  ////hctau3D_eff_1->GetYaxis()->CenterTitle();
  //hctau3D_eff_1->GetYaxis()->SetTitle("Efficiency");
  //hctau3D_eff_1->GetYaxis()->SetRangeUser(0.,1.1);
  ////hctau3D_eff_1->Draw("E");
  ////eeeeee
  //hctau3D_eff_2->SetMarkerStyle(24);
  //hctau3D_eff_2->SetMarkerColor(1);
  //hctau3D_eff_2->SetLineColor(1);
  //hctau3D_eff_2->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  //hctau3D_eff_2->GetYaxis()->SetTitle("Efficiency");
  //hctau3D_eff_2->GetYaxis()->SetRangeUser(0.,1.3);
  //hctau3D_eff_2->Draw("E");
  //lt1->SetTextSize(0.03);
  //if(state==1) lt1->DrawLatex(0.13,0.84,"Prompt J/#psi (PbPb)");
  //else if(state==2) lt1->DrawLatex(0.13,0.84,"Non-Prompt J/#psi (PbPb)");
  //lt1->DrawLatex(0.13,0.80,"|y| < 2.4");
  //lt1->SetTextSize(0.02);
  ////lt1->DrawLatex(0.15,0.88,Form("Yields : %0.2f #pm %0.2f",yields[0], yields[1]));
  ////hctau3D_eff_2->SetMarkerStyle(25);
  ////hctau3D_eff_2->SetMarkerColor(2);
  ////hctau3D_eff_2->SetLineColor(2);
  ////hctau3D_eff_2->Draw("same");
  //hctau3D_eff_3->SetMarkerStyle(30);
  //hctau3D_eff_3->SetMarkerColor(kOrange+1);
  //hctau3D_eff_3->SetLineColor(kOrange+1);
  //hctau3D_eff_3->Draw("same");
  //hctau3D_eff_4->SetMarkerStyle(26);
  //hctau3D_eff_4->SetMarkerColor(3);
  //hctau3D_eff_4->SetLineColor(3);
  //hctau3D_eff_4->Draw("same");
  //hctau3D_eff_5->SetMarkerStyle(27);
  //hctau3D_eff_5->SetMarkerColor(4);
  //hctau3D_eff_5->SetLineColor(4);
  //hctau3D_eff_5->Draw("same");
  ////cctau3D_eff->BuildLegend();
  ////legend->AddEntry("hctau3D_eff_1","|y|: 0-2.4, 0-100%","lep");
  ////legend->AddEntry("hctau3D_eff_2","|y|: 0-0.6, 0-100%","lep");
  ////legend->AddEntry("hctau3D_eff_3","|y|: 0.6-1.2, 0-100%","lep");
  ////legend->AddEntry("hctau3D_eff_4","|y|: 1.2-1.8, 0-100%","lep");
  ////legend->AddEntry("hctau3D_eff_5","|y|: 1.8-2.4, 0-100%","lep");
  ////legend->SetBorderSize( 0);
  ////legend->AddEntry("hctau3D_eff_1","0-90%","lep");
  //legend->AddEntry("hctau3D_eff_2","0-10%","lep");
  //legend->AddEntry("hctau3D_eff_3","10-30%","lep");
  //legend->AddEntry("hctau3D_eff_4","30-60%","lep");
  //legend->AddEntry("hctau3D_eff_5","60-90%","lep");
  //legend->SetBorderSize(0);
  ////legend->SetFillColor ( 0);
  //legend->Draw("same");
  ////legend->Draw();
  ////cctau3D_eff->SaveAs("Eff_pt.png");


  //Save efficiency files for later use.
  hctau3D_eff_1->SetName(Form("mc_eff_vs_pt_TnP%d_1 ",isTnP));
  hctau3D_eff_2->SetName(Form("mc_eff_vs_pt_TnP%d_2  ",isTnP));
  hctau3D_eff_3->SetName(Form("mc_eff_vs_pt_TnP%d_3 ",isTnP));
  hctau3D_eff_4->SetName(Form("mc_eff_vs_pt_TnP%d_4",isTnP));
  hctau3D_eff_5->SetName(Form("mc_eff_vs_pt_TnP%d_5",isTnP));
  //TString outFileName = Form("mc_eff_vs_pt_TnP%d_PtW%d_OfficialMC_Y%dS_muPtCut%.1f.root",isTnP,isPtWeight,state,muPtCut);
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hctau3D_eff_1->Write();
  hctau3D_eff_2->Write();
  hctau3D_eff_3->Write();
  hctau3D_eff_4->Write();
  hctau3D_eff_5->Write();
  hctau3D_gen_1->Write();
  hctau3D_gen_2->Write();
  hctau3D_gen_3->Write();
  hctau3D_gen_4->Write();
  hctau3D_gen_5->Write();
  hctau3D_reco_1->Write();
  hctau3D_reco_2->Write();
  hctau3D_reco_3->Write();
  hctau3D_reco_4->Write();
  hctau3D_reco_5->Write();
  //if(isTnP) hctau3D_tnp_trig->Write();
  outFile->Close();

}
