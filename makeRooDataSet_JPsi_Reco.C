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

void makeRooDataSet_JPsi_Reco(bool isMC = true, bool fPtW = true, bool fAccW = true, bool fEffW = true, int state=1) //state 0: inclusive, state 1: Prompt, state 2: NonPrompt
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
    if(state==1) rf = new TFile("skimmedFiles/OniaFlowSkim_Jpsi_MC_Prompt_210107.root","read");
    else if(state==2) rf = new TFile("skimmedFiles/OniaFlowSkim_BtoJpsi_GENONLY_210128.root","read");
  }
  //else if(!isMC) rf = new TFile("./skimmedFiles/OniaFlowSkim_JpsiTrig_AllPD_isMC0_HFNom_200407.root","read");
  else if(!isMC) rf = new TFile("/work2/Oniatree/JPsi/skimmed_file/OniaFlowSkim_BtoJpsi_GENONLY_210128.root","read");
  //  else if(!isMC) rf = new TFile("OniaFlowSkim_JpsiTrig_JPsi_isMC0_HFNom_200619_RooDataSet_test100k.root","read");
  TTree *tree = (TTree*) rf -> Get("mmepevt");
  
  TFile *fPtW1 = new TFile("Eff_Acc/WeightedFcN_fit/ratioDataMC_AA_PromptJpsi_DATA_All_y.root","read");
  TFile *fPtW2 = new TFile("Eff_Acc/WeightedFcN_fit/ratioDataMC_AA_PromptJpsi_DATA_Forward_y.root","read");
  TF1* fptw1 = (TF1*) fPtW1->Get("dataMC_Ratio1");
  TF1* fptw2 = (TF1*) fPtW2->Get("dataMC_Ratio1");

  bool isTnP = true;
  bool isPtW = true;
  TFile *fEff = new TFile(Form("./Eff_Acc/mc_eff_vs_pt_cent_20_to_120_rap_%s_pbpb_Jpsi_PtW%d_tnp%d_drawsame2.root",bCont.Data(),isPtW,isTnP),"read");
  TH1D* hEffPt[15];
  hEffPt[0] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy0_1p2",isTnP,isPtW));
  hEffPt[1] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p2_1p6",isTnP,isPtW));
  hEffPt[2] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p6_1p8",isTnP,isPtW));
  hEffPt[3] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_20_to_120_absy1p8_2p4",isTnP,isPtW));
  
  TFile *fAcc = new TFile("./Eff_Acc/Acc_WeightVariation_20210426.root","read");
  TH1D* hAccPt[3];
  hAccPt[0] = (TH1D*) fAcc -> Get("hAccPt_2021_ally_wt");
  hAccPt[1] = (TH1D*) fAcc -> Get("hAccPt_2021_midy_wt");
  hAccPt[2] = (TH1D*) fAcc -> Get("hAccPt_2021_Fory_wt");

  //SetBranchAddress
  const int nMaxDimu = 1000;
  float mass[nMaxDimu];
  float pt[nMaxDimu];
  float y[nMaxDimu]; 
  float pt1[nMaxDimu];
  float pt2[nMaxDimu];
  float eta[nMaxDimu];
  float eta1[nMaxDimu];
  float eta2[nMaxDimu];
  float ctau3D[nMaxDimu];
  Int_t event; 
  Int_t nDimu; 
  double weight;

  TBranch *b_event;
  TBranch *b_nDimu;
  TBranch *b_mass;
  TBranch *b_pt;
  TBranch *b_y;
  TBranch *b_eta;
  TBranch *b_eta1;
  TBranch *b_eta2;
  TBranch *b_ctau3D;
  TBranch *b_pt1;
  TBranch *b_pt2;
  TBranch *b_weight;

  tree -> SetBranchAddress("event", &event, &b_event);
  tree -> SetBranchAddress("nDimu", &nDimu, &b_nDimu);
  tree -> SetBranchAddress("mass", mass, &b_mass);
  tree -> SetBranchAddress("y", y, &b_y);
  tree -> SetBranchAddress("pt", pt, &b_pt);
  tree -> SetBranchAddress("pt1", pt1, &b_pt1);
  tree -> SetBranchAddress("pt2", pt2, &b_pt2);
  tree -> SetBranchAddress("eta", eta, &b_eta);
  tree -> SetBranchAddress("eta1", eta1, &b_eta1);
  tree -> SetBranchAddress("eta2", eta2, &b_eta2);
  tree -> SetBranchAddress("ctau3D", ctau3D, &b_ctau3D);
  tree -> SetBranchAddress("weight", &weight, &b_weight);

  ////////////////////////////////////////////////////////////////////////
  //////////////////  RooDataSet 
  ////////////////////////////////////////////////////////////////////////
  RooRealVar* massVar  = new RooRealVar("mass","mass variable",1.0,6.0,"GeV/c^{2}");
  RooRealVar* ptVar    = new RooRealVar("pt","pt variable", 0,100,"GeV/c");
  RooRealVar* yVar     = new RooRealVar("y","rapidity of the dimuon pair", -5,5,"");
  RooRealVar* pt1Var   = new RooRealVar("pt1","pt of muon+", 0,500,"GeV/c");
  RooRealVar* eta1Var  = new RooRealVar("eta1","eta of muon+", -4,4,"");
  RooRealVar* pt2Var   = (RooRealVar*)pt1Var->Clone("pt2");
  RooRealVar* eta2Var  = (RooRealVar*)eta1Var->Clone("eta2");
  RooRealVar* evtWeight = new RooRealVar("weight","corr weight", 0, 10000,"");
  RooRealVar* ctau3DVar = new RooRealVar("ctau3D","ctau3D dimuon",-100000.0, 100000.0,"mm");
  RooRealVar* NumDimu = new RooRealVar("NumDimu","number of dimuon",0,100,"");
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var,*evtWeight);
  argSet->add(*NumDimu); argSet->add(*ctau3DVar);

  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);

  int nDimuPass=0;
  int nDimu_one=0;
  int nDimu_more=0;
  int nDimu_all=0;

  double weight_pt = 1;
  double weight_acc = 1;
  double weight_eff = 1;
  double weight_;

  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;

  //  double SiMuPtCut = 3.5;
  //Begin Loop
  for(int i=0; i<nEvt; i++){
    tree->GetEntry(i);
    nDimuPass=0;

    nDimu_all++;
    if(nDimuPass>1) {nDimu_more++; continue;}
    if(nDimuPass==1) nDimu_one++;

    // Fill Dimuon Loop
    for(int j=0; j<nDimu; j++){
      if(pt[j]>=3 && pt[j]<50){
        //if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && (double)pt1[j]>SiMuPtCut&&(double)pt2[j]>SiMuPtCut&&abs((double)eta1[j])<2.4&&abs((double)eta2[j])<2.4)) continue;
        if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && IsAcceptanceQQ(pt1[j],eta1[j])&&IsAcceptanceQQ(pt2[j],eta2[j])) ) continue;

        if(fPtW){
          if(fabs(abs((double)y[j])) < 1.6 ) weight_pt = fptw1->Eval(pt[j]);
          else if(fabs(abs((double)y[j])) > 1.6 ) weight_pt = fptw2->Eval(pt[j]);
        }
        if(fAccW){
          if ( abs((double)y[j])<1.6) { weight_acc = getAccWeight(hAccPt[1], pt[j]); }
          else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<2.4 ) { weight_acc = getAccWeight(hAccPt[2], pt[j]); }
        }
        if(fEffW){
          if ( abs((double)y[j])<1.2) { weight_eff = getEffWeight(hEffPt[0], pt[j]); }
          else if ( abs((double)y[j])>=1.2 && abs((double)y[j])<1.6 ) { weight_eff = getEffWeight(hEffPt[1], pt[j]); }
          else if ( abs((double)y[j])>=1.6 && abs((double)y[j])<1.8 ) { weight_eff = getEffWeight(hEffPt[2], pt[j]); }
          else if ( abs((double)y[j])>=1.8 && abs((double)y[j])<2.4 ) { weight_eff = getEffWeight(hEffPt[3], pt[j]); }
        }

        weight_ = weight*weight_pt*weight_acc*weight_eff;
        //cout << weight << " " << weight_eff << " " << weight_acc << " " << weight_ << endl;
        massVar->setVal( (double)mass[j] ) ;
        ptVar->setVal(   (double)pt[j]   ) ;
        yVar->setVal(    (double)y[j]    ) ;
        pt1Var->setVal(  (double)pt1[j]  ) ;
        eta1Var->setVal( (double)eta1[j] ) ;
        pt2Var->setVal(  (double)pt2[j]  ) ;
        eta2Var->setVal( (double)eta2[j] ) ;
        ctau3DVar->setVal( (double)ctau3D[j] ) ;
        evtWeight->setVal( (double)weight_ ) ;
        NumDimu->setVal((int)nDimu);
        dataSet->add( *argSet);
      }
    }
  }

  cout << "All : " << nDimu_all << endl;
  cout << "more than one dimuon : " << nDimu_more << endl;
  cout << "one dimuon : " << nDimu_one << endl;

  if (isMC && state==1) {TFile *wf = new TFile(Form("skimmedFiles/OniaRooDataSet_isMC%d_PR_JPsi_20210601.root",isMC),"recreate");  wf->cd();}
  else if (isMC && state==2) {TFile *wf = new TFile("skimmedFiles/OniaRooDataSet_NonPrompt_GenInReco.root","recreate");  wf->cd();}

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
