#include <iostream>
#include "commonUtility.h"
#include "cutsAndBin.h"
#include "HiEvtPlaneList.h"
#include "Style.h"
#include "tdrstyle.C"
#include "CMS_lumi_v2mass.C"
#include "rootFitHeaders.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"

using namespace std;
using namespace RooFit;

using namespace hi;

double getAccWeight(TH1D* h = 0, double pt = 0);
double getEffWeight(TH1D* h = 0, double pt = 0);
void GetHistSqrt(TH1D* h1 =0, TH1D* h2=0);

void makeRooDataSet_JPsi(bool isMC = false, bool fAccW = true, bool fEffW = true, int state=2) //state 0: inclusive, state 1: Prompt, state 2: NonPrompt
{
  //Basic Setting
  gStyle->SetOptStat(0);
  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  TString dimusignString;
  TString bCont;
  TString outName;
  if (state==1) { bCont = "prompt"; outName="PR"; }
  else if (state==2) {bCont = "nprompt"; outName="NP"; }

  //READ Input Skimmed File
  TFile *rf;
  if(isMC){
    if(state==1) rf = new TFile("/work2/Oniatree/JPsi/skimmed_file/OniaFlowSkim_JpsiTrig_isMC1_Prompt_HFNom_200407.root","read");
    else if(state==2) rf = new TFile("/work2/Oniatree/JPsi/skimmed_file/OniaFlowSkim_JpsiTrig_isMC1_NonPrompt_HFNom_200407.root","read");
  }
  //else if(!isMC) rf = new TFile("./skimmedFiles/OniaFlowSkim_JpsiTrig_AllPD_isMC0_HFNom_200407.root","read");
  else if(!isMC) rf = new TFile("./skimmedFiles/OniaFlowSkim_JpsiTrig_DBAllPD_isMC0_HFNom_201127.root","read");
//  else if(!isMC) rf = new TFile("OniaFlowSkim_JpsiTrig_JPsi_isMC0_HFNom_200619_RooDataSet_test100k.root","read");
  TTree *tree = (TTree*) rf -> Get("mmepevt");

  
  //Get Correction histograms
  bool isTnP = true;
  bool isPtW = true;
//  TFile *fEff = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/Efficiency/mc_eff_vs_pt_TnP%d_PtW1_OfficialMC_Y%dS_muPtCut3.5.root",isTnP,state),"read");
  TFile *fEff = new TFile(Form("mc_eff_vs_pt_cent_rap_%s_pbpb_Jpsi.root",bCont.Data()),"read");
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
  TFile *fAcc = new TFile(Form("mc_Acc_vs_pt_cent_rap_%s_pbpb_Jpsi_PtW1.root",bCont.Data()),"read");
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
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  RooDataSet 
  ////////////////////////////////////////////////////////////////////////
  RooRealVar* massVar  = new RooRealVar("mass","mass variable",0,200,"GeV/c^{2}");
  RooRealVar* ptVar    = new RooRealVar("pt","pt variable", 0,100,"GeV/c");
  RooRealVar* yVar     = new RooRealVar("y","rapidity of the dimuon pair", -5,5,"");
  RooRealVar* pt1Var   = new RooRealVar("pt1","pt of muon+", 0,500,"GeV/c");
  RooRealVar* eta1Var  = new RooRealVar("eta1","eta of muon+", -4,4,"");
  RooRealVar* pt2Var   = (RooRealVar*)pt1Var->Clone("pt2");
  RooRealVar* eta2Var  = (RooRealVar*)eta1Var->Clone("eta2");
  RooRealVar* cBinVar   = new RooRealVar("cBin","Centrality bin", -100,500,"");
  RooRealVar* ep2Var   = new RooRealVar("ep2","2nd order event plane", -100,100,"");
  RooRealVar* evtWeight = new RooRealVar("weight","corr weight", 0, 10000,"");
  RooRealVar* recoQQ = new RooRealVar("recoQQsign","qq sign",-1,3,"");
  RooRealVar* ctau3DVar = new RooRealVar("ctau3D","ctau3D dimuon",-20,20,"cm");
  RooRealVar* ctau3DErrVar = new RooRealVar("ctau3DErr","ctau3DErr dimuon",0,0.5,"cm");
  RooRealVar* ctau3DResVar = new RooRealVar("ctau3DRes","ctau3D Resolution dimuon",-20,20,"cm");
  RooRealVar* ctau3DSignVar = new RooRealVar("ctau3DSign","ctau3D Significance dimuon",-20,20,"cm");
  RooRealVar* NumDimu = new RooRealVar("NumDimu","number of dimuon",0,100,"");
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var,*evtWeight);
  argSet->add(*cBinVar); argSet->add(*recoQQ); argSet->add(*NumDimu); argSet->add(*ctau3DVar); argSet->add(*ctau3DErrVar); argSet->add(*ctau3DResVar);
  
  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);


  int nDimuPass=0;
  int nDimu_one=0;
  int nDimu_more=0;
  int nDimu_all=0;

  double weight_acc = 1;
  double weight_eff = 1;
  

  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;
  
  double SiMuPtCut = 3.5;
  //Begin Loop
  for(int i=0; i<nEvt; i++){
    tree->GetEntry(i);
    nDimuPass=0;
    
    if(fabs(vz)>=15) continue;

	
    //Remove double candidate
    for(int j=0; j<nDimu; j++){
      //if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && (double)pt1[j]>SiMuPtCut&&(double)pt2[j]>SiMuPtCut&&abs((double)eta1[j])<2.4&&abs((double)eta2[j])<2.4)) continue;
      if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && IsAcceptanceQQ(pt1[j],eta1[j]) && IsAcceptanceQQ(pt2[j],eta2[j])) ) continue;
      nDimuPass++;
    }


    nDimu_all++;
    if(nDimuPass>1) {nDimu_more++; continue;}
    if(nDimuPass==1) nDimu_one++;

    // Fill Dimuon Loop
    for(int j=0; j<nDimu; j++){
        //if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && (double)pt1[j]>SiMuPtCut&&(double)pt2[j]>SiMuPtCut&&abs((double)eta1[j])<2.4&&abs((double)eta2[j])<2.4)) continue;
        if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && IsAcceptanceQQ(pt1[j],eta1[j])&&IsAcceptanceQQ(pt2[j],eta2[j])) ) continue;
        weight_acc = 1;
        weight_eff = 1;
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
        recoQQ->setVal((int)recoQQsign[j]);     
        massVar->setVal( (double)mass[j] ) ;
        ptVar->setVal(   (double)pt[j]   ) ;
        yVar->setVal(    (double)y[j]    ) ;
        pt1Var->setVal(  (double)pt1[j]  ) ;
        eta1Var->setVal( (double)eta1[j] ) ;
        pt2Var->setVal(  (double)pt2[j]  ) ;
        eta2Var->setVal( (double)eta2[j] ) ;
        cBinVar->setVal( (double)cBin ) ;
		ctau3DVar->setVal( (double)ctau3D[j] ) ;
		ctau3DErrVar->setVal( (double)ctau3DErr[j] ) ;
		ctau3DResVar->setVal( (double)ctau3D[j]/ctau3DErr[j] ) ;
		ctau3DSignVar->setVal( (double)ctau3DErr[j]/ctau3D[j] ) ;
        evtWeight->setVal( (double)weight_ ) ;
        NumDimu->setVal((int)nDimu);
        dataSet->add( *argSet);
    }
  }

  cout << "All : " << nDimu_all << endl;
  cout << "more than one dimuon : " << nDimu_more << endl;
  cout << "one dimuon : " << nDimu_one << endl;
  
  
  if (isMC && state==1) {TFile *wf = new TFile(Form("skimmedFiles/OniaRooDataSet_isMC%d_PR_JPsi_20201123.root",isMC),"recreate");  wf->cd();}
  else if (isMC && state==2) {TFile *wf = new TFile(Form("skimmedFiles/OniaRooDataSet_isMC%d_BtoJPsi_20201123.root",isMC),"recreate");  wf->cd();}
  else if (!isMC) {TFile *wf = new TFile(Form("skimmedFiles/OniaRooDataSet_isMC%d_JPsi_%sw_Effw%d_Accw%d_PtW%d_TnP%d_20201216.root",isMC,outName.Data(),fEffW,fAccW,isPtW,isTnP),"recreate");  wf->cd();}
 dataSet->Write();
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
