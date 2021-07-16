#include <ctime>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBin.h"
//#include "tnp_weight_lowptPbPb.h"
#include "Eff_Acc/tnp_weight_lowptPbPb_num_den.h"

static const long MAXTREESIZE = 1000000000000;
double getAccWeight(TH1D* h = 0, double pt = 0);
double getEffWeight(TH1D* h = 0, double pt = 0);


void SkimTree_Event_Jpsi_GH_test_v3(int nevt=-1, bool isMC = true, int kTrigSel = kTrigJpsi, int hiHFBinEdge = 0, int PDtype = 1) 
{

  using namespace std;
  using namespace hi;

  // Example of using event plane namespace 
  cout << " Index of "<< EPNames[HFm2] << " = " << HFm2 << endl;
  cout << " Index of "<< EPNames[HFp2] << " = " << HFp2 << endl;
  cout << " Index of "<< EPNames[trackmid2] << " = " << trackmid2 << endl;

  TString fnameData1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v1_Oniatree_addvn_part*.root";
  TString fnameData2 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v2_Oniatree_addvn_part*.root";
  TString fnameDataReReco = "/work2/Oniatree/ReReco/DoubleMuon/ReReco_Oniatree_addvn_part*.root";
  TString fnameDataReRecoPeri = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/ReReco/AOD/DoubleMuonPsiPeri/ReReco_Oniatree_addvn_part*.root";
  //TString fnameMC = "/eos/cms/store/group/phys_heavyions/gbak/2018PbPbMC/JPsi/OniatreeMC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8.root";
  //TString fnameMC = "/eos/cms/store/group/phys_heavyions/gbak/2018PbPbMC/JPsi/OniatreeMC_JPsi_pThat2_TunedCP5_HydjetDrumMB_5p02TeV_Pythia8_part*.root";
  TString fnameMC = "/work2/Oniatree/JPsi/OniatreeMC_JPsi_pThat2_TunedCP5_HydjetDrumMB_5p02TeV_Pythia8_part*.root";

  TString fPD;
  if(PDtype==1) fPD = "DB";
  else if(PDtype==2) fPD = "DBPeri";

  TChain *mytree = new TChain("myTree");
  if(!isMC){
    if(PDtype==1) mytree->Add(fnameDataReReco.Data());
    else if(PDtype==2) mytree->Add(fnameDataReRecoPeri.Data());
  }
  else if(isMC){
    mytree->Add(fnameMC.Data());
  }

  const int maxBranchSize = 1000;

  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Float_t         SumET_HF;
  Int_t           Reco_QQ_size;
  Int_t           Reco_mu_size;
  Int_t		  Gen_QQ_size;
  Int_t 	  Gen_mu_size;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_mu_4mom;
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_mu_4mom;
  ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
  ULong64_t       Reco_mu_trig[maxBranchSize];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_runNb;   //!
  TBranch        *b_eventNb;   //!
  TBranch        *b_LS;
  TBranch        *b_zVtx;   //!
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_SumET_HF;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_mu_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_mu_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_mu_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!
  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_mu_size;   //!
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Gen_mu_4mom;   //!

  Bool_t          Reco_mu_highPurity[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_highPurity;   //!
  mytree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);

  Reco_QQ_4mom = 0;
  Reco_mu_4mom = 0;
  Gen_QQ_4mom = 0;
  Gen_mu_4mom = 0;
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  mytree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  mytree->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
  mytree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  mytree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  //  muon id 
  Int_t           Reco_QQ_mupl_idx[maxBranchSize];
  Int_t           Reco_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Reco_QQ_mupl_idx;
  TBranch        *b_Reco_QQ_mumi_idx;
  mytree->SetBranchAddress("Reco_QQ_mupl_idx",Reco_QQ_mupl_idx,&b_Reco_QQ_mupl_idx);
  mytree->SetBranchAddress("Reco_QQ_mumi_idx",Reco_QQ_mumi_idx,&b_Reco_QQ_mumi_idx);

  Int_t           Gen_QQ_mupl_idx[maxBranchSize];
  Int_t           Gen_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Gen_QQ_mupl_idx;
  TBranch        *b_Gen_QQ_mumi_idx;
  mytree->SetBranchAddress("Gen_QQ_mupl_idx",Gen_QQ_mupl_idx,&b_Gen_QQ_mupl_idx);
  mytree->SetBranchAddress("Gen_QQ_mumi_idx",Gen_QQ_mumi_idx,&b_Gen_QQ_mumi_idx);

  Int_t           Reco_mu_nTrkHits[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkHits;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  Float_t         Reco_mu_normChi2_global[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_normChi2_global;   //!
  mytree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  Int_t           Reco_mu_nMuValHits[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nMuValHits;   //!
  mytree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  Int_t           Reco_mu_StationsMatched[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_StationsMatched;   //!
  mytree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  Float_t         Reco_mu_dxy[maxBranchSize];   //[Reco_mu_size]
  Float_t         Reco_mu_dxyErr[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dxy;   //!
  TBranch        *b_Reco_mu_dxyErr;   //!
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  mytree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  Float_t         Reco_mu_dz[maxBranchSize];   //[Reco_mu_size]
  Float_t         Reco_mu_dzErr[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  Int_t           Reco_mu_nTrkWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Bool_t          Reco_mu_TMOneStaTight[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_TMOneStaTight;   //!

  mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  Int_t           Reco_mu_nPixWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Int_t           Reco_QQ_sign[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  Float_t         rpAng[29];   //[nEP]
  TBranch        *b_rpAng;   //!
  //  mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

  Int_t           Reco_mu_nPixValHits[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  Float_t         Reco_mu_ptErr_global[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_ptErr_global;   //!
  mytree->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);

  Int_t           Reco_mu_SelectionType[maxBranchSize];
  TBranch        *b_Reco_mu_SelectionType;
  mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

  Float_t         Reco_QQ_ctau3D[maxBranchSize];
  Float_t         Reco_QQ_ctauErr3D[maxBranchSize];
  TBranch        *b_Reco_QQ_ctau3D;
  TBranch        *b_Reco_QQ_ctauErr3D;
  mytree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
  mytree->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);

  Float_t         Gen_QQ_ctau3D[maxBranchSize];
  TBranch        *b_Gen_QQ_ctau3D;
  mytree->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);

  Int_t Reco_QQ_whichGen[maxBranchSize];
  Int_t Gen_QQ_whichRec[maxBranchSize];
  Int_t Reco_mu_whichGen[maxBranchSize];
  Int_t Gen_mu_whichRec[maxBranchSize];
  TBranch *b_Reco_QQ_whichGen;
  TBranch *b_Gen_QQ_whichRec;
  TBranch *b_Reco_mu_whichGen;
  TBranch *b_Gen_mu_whichRec;
  Float_t Gen_weight;
  TBranch *b_Gen_weight;
  if(isMC){
    mytree->SetBranchAddress("Reco_QQ_whichGen",Reco_QQ_whichGen, &b_Reco_QQ_whichGen);
    mytree->SetBranchAddress("Gen_QQ_whichRec",Gen_QQ_whichRec, &b_Gen_QQ_whichRec);
    mytree->SetBranchAddress("Reco_mu_whichGen",Reco_mu_whichGen, &b_Reco_mu_whichGen);
    mytree->SetBranchAddress("Gen_mu_whichRec",Gen_mu_whichRec, &b_Gen_mu_whichRec);
    mytree->SetBranchAddress("Gen_weight",&Gen_weight, &b_Gen_weight);
  }

  TChain *eptree = new TChain("tree");
  if(!isMC){
    eptree->Add(fnameDataReReco.Data());
  }
  else if(isMC){
    eptree->Add(fnameMC.Data());
  }

  const int nEP = 29;  // number of event planes in the tree
  double qx[nEP]; 
  double qy[nEP]; 
  TBranch *b_qx;
  TBranch *b_qy;
  eptree->SetBranchAddress("qx",qx, &b_qx);
  eptree->SetBranchAddress("qy",qy, &b_qy);

  int trigIndx=0;
  if(kTrigSel == kTrigJpsi) trigIndx=0;
  else if(kTrigSel == kTrigUps) trigIndx=1;
  else if(kTrigSel == kTrigL1DBOS40100) trigIndx=2;
  else if(kTrigSel == kTrigL1DB50100) trigIndx=3;

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

  int kL2filter = 16;
  int kL3filter = 17;

  int count =0;
  int counttnp =0;

  TString fCentSelHF = "HFNom";
  if(hiHFBinEdge==1) fCentSelHF = "HFUp";
  else if(hiHFBinEdge==-1) fCentSelHF = "HFDown";
  TFile* newfile;
  //newfile = new TFile(Form("OniaFlowSkim_%sTrig_%sPD_isMC%d_%s_200109.root",fTrigName[trigIndx].Data(),fPD.Data(),isMC,fCentSelHF.Data()),"recreate");
  newfile = new TFile(Form("OniaFlowSkim_%sTrig_JPsi_isMC%d_%s_CtauResSys_210714.root",fTrigName[trigIndx].Data(),isMC,fCentSelHF.Data()),"recreate");

  const static int nMaxDimu = 1000;
  int evt;
  int runN;
  int lumi;
  int cBin;
  int nDimu;
  int nDimuGen;
  float vz;
  int recoQQsign[nMaxDimu];
  float mass[nMaxDimu];
  float pt[nMaxDimu];
  float pt1[nMaxDimu];
  float pt2[nMaxDimu];
  float y[nMaxDimu];
  float eta[nMaxDimu];
  float eta1[nMaxDimu];
  float eta2[nMaxDimu];
  float phi[nMaxDimu];
  float phi1[nMaxDimu];
  float phi2[nMaxDimu];
  float qxa[nMaxDimu];
  float qya[nMaxDimu];
  float qxb[nMaxDimu];
  float qyb[nMaxDimu];
  float qxc[nMaxDimu];
  float qyc[nMaxDimu];
  float qxdimu[nMaxDimu];
  float qydimu[nMaxDimu];
  float qxmupl[nMaxDimu];
  float qxmumi[nMaxDimu];
  float qymupl[nMaxDimu];
  float qymumi[nMaxDimu];
  float ctau3D[nMaxDimu];
  float ctau3DErr[nMaxDimu];
  float ctau3Dtrue[nMaxDimu];
  double TnPweight[nMaxDimu] = {1.};
  float Genmass[nMaxDimu];
  float Genpt[nMaxDimu];
  float Genpt1[nMaxDimu];
  float Genpt2[nMaxDimu];
  float Geny[nMaxDimu];
  float Geneta[nMaxDimu];
  float Geneta1[nMaxDimu];
  float Geneta2[nMaxDimu];
  float Genphi[nMaxDimu];
  float Genphi1[nMaxDimu];
  float Genphi2[nMaxDimu];
  float Genctau3D[nMaxDimu];
  double weight = 1;

  TTree* mmevttree = new TTree("mmepevt","dimuonAndEventPlanes in event based");
  mmevttree->SetMaxTreeSize(MAXTREESIZE);
  mmevttree->Branch("event",&evt,"event/I");
  mmevttree->Branch("runN",&runN,"runN/I");
  mmevttree->Branch("lumi",&lumi,"lumi/I");
  mmevttree->Branch("cBin",&cBin,"cBin/I");
  mmevttree->Branch("vz",&vz,"vz/F");
  mmevttree->Branch("nDimu",&nDimu,"nDimu/I");
  mmevttree->Branch("mass",mass,"mass[nDimu]/F");
  mmevttree->Branch("y",y,"y[nDimu]/F");
  mmevttree->Branch("pt",pt,"pt[nDimu]/F");
  mmevttree->Branch("pt1",pt1,"pt1[nDimu]/F");
  mmevttree->Branch("pt2",pt2,"pt2[nDimu]/F");
  mmevttree->Branch("eta",eta,"eta[nDimu]/F");
  mmevttree->Branch("eta1",eta1,"eta1[nDimu]/F");
  mmevttree->Branch("eta2",eta2,"eta2[nDimu]/F");
  mmevttree->Branch("phi",phi,"phi[nDimu]/F");
  mmevttree->Branch("phi1",phi1,"phi1[nDimu]/F");
  mmevttree->Branch("phi2",phi2,"phi2[nDimu]/F");
  mmevttree->Branch("qxa",qxa,"qxa[nDimu]/F");
  mmevttree->Branch("qxb",qxb,"qxb[nDimu]/F");
  mmevttree->Branch("qxc",qxc,"qxc[nDimu]/F");
  mmevttree->Branch("qya",qya,"qya[nDimu]/F");
  mmevttree->Branch("qyb",qyb,"qyb[nDimu]/F");
  mmevttree->Branch("qyc",qyc,"qyc[nDimu]/F");
  mmevttree->Branch("qxdimu",qxdimu,"qxdimu[nDimu]/F");
  mmevttree->Branch("qydimu",qydimu,"qydimu[nDimu]/F");
  mmevttree->Branch("qxmupl",qxmupl,"qxmupl[nDimu]/F");
  mmevttree->Branch("qxmumi",qxmumi,"qxmumi[nDimu]/F");
  mmevttree->Branch("qymupl",qymupl,"qymupl[nDimu]/F");
  mmevttree->Branch("qymumi",qymumi,"qymumi[nDimu]/F");
  mmevttree->Branch("recoQQsign",recoQQsign,"recoQQsign[nDimu]/I");
  mmevttree->Branch("ctau3D",ctau3D,"ctau3D[nDimu]/F");
  mmevttree->Branch("ctau3DErr",ctau3DErr,"ctau3DErr[nDimu]/F");
  mmevttree->Branch("ctau3Dtrue",ctau3Dtrue,"ctau3Dtrue[nDimu]/F");
  mmevttree->Branch("TnPweight",TnPweight,"TnPweight[nDimu]/D");
  mmevttree->Branch("weight",&weight,"weight/D");

  TTree* mmgentree = new TTree("mmGentree","Gen Di-muon Pairs");
  mmgentree->SetMaxTreeSize(MAXTREESIZE);
  mmgentree->Branch("event",&evt,"event/I");
  mmgentree->Branch("nDimuGen",&nDimuGen,"nDimuGen/I");
  mmgentree->Branch("mass",Genmass,"Genmass[nDimuGen]/F");
  mmgentree->Branch("y",Geny,"Geny[nDimuGen]/F");
  mmgentree->Branch("pt",Genpt,"Genpt[nDimuGen]/F");
  mmgentree->Branch("pt1",Genpt1,"Genpt1[nDimuGen]/F");
  mmgentree->Branch("pt2",Genpt2,"Genpt2[nDimuGen]/F");
  mmgentree->Branch("eta",Geneta,"Geneta[nDimuGen]/F");
  mmgentree->Branch("eta1",Geneta1,"Geneta1[nDimuGen]/F");
  mmgentree->Branch("eta2",Geneta2,"Geneta2[nDimuGen]/F");
  mmgentree->Branch("phi",Genphi,"Genphi[nDimuGen]/F");
  mmgentree->Branch("phi1",Genphi1,"Genphi1[nDimuGen]/F");
  mmgentree->Branch("phi2",Genphi2,"Genphi2[nDimuGen]/F");
  mmgentree->Branch("ctau3D",Genctau3D,"Genctau3D[nDimuGen]/F");

  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;

  TLorentzVector* JP_Gen = new TLorentzVector;
  TLorentzVector* mupl_Gen = new TLorentzVector;
  TLorentzVector* mumi_Gen = new TLorentzVector;
  ////////////////////////////////////////////////////////////////////////
  ////////////////// RooDataSet 
  ////////////////////////////////////////////////////////////////////////
  //RooRealVar* Ctau3DErrVar = new RooRealVar("Ctau3DErr","Ctau Error variable",0,350,"");
  //RooArgSet* argSet = new RooArgSet(*Ctau3DErr);

  //RooDataSet* dataSet = new RooDataSet("dataSet", "a dataset", argSet);


  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();

  cout << "Total events = " << nevt << ", : " << eptree->GetEntries() << endl;
  //////////////////////////////////////////////////////
  //////////////////////// Reco ////////////////////////
  //////////////////////////////////////////////////////

  Int_t mupl_idx;
  Int_t mumi_idx;
  cout<< " " << endl;
  cout<< "Start fill Reco tree"<<endl;
  for(int iev=0; iev<nevt ; ++iev)
  {
    //if(iev==10000)break;
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
    eptree->GetEntry(iev);

    nDimu = 0;

    if(PDtype == 1 && runNb >= 327123) continue;

    if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      runN = runNb;
      evt = eventNb;
      lumi = LS;
      cBin = -999;
      if(hiHFBinEdge ==0) cBin = getHiBinFromhiHF(SumET_HF);
      else if(hiHFBinEdge == 1) cBin = getHiBinFromhiHF_Up(SumET_HF);
      else if(hiHFBinEdge == -1) cBin = getHiBinFromhiHF_Down(SumET_HF);
      if(cBin==-999){ cout << "ERROR!!! No HF Centrality Matching!!" << endl; return;}
      vz = zVtx;

      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

      weight = 1.;
      if(isMC) weight = findNcoll(Centrality) * Gen_weight;

      if(isMC){
        if(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
        if(Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;
        mupl_Gen = (TLorentzVector*) Gen_mu_4mom->At(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]]);
        mumi_Gen = (TLorentzVector*) Gen_mu_4mom->At(Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]]);
        if(sqrt(pow(mupl_Reco->Eta()-mupl_Gen->Eta(),2)+pow(mupl_Reco->Phi()-mupl_Gen->Phi(),2)>0.03)) continue;
        if(sqrt(pow(mumi_Reco->Eta()-mumi_Gen->Eta(),2)+pow(mumi_Reco->Phi()-mumi_Gen->Phi(),2)>0.03)) continue;
      }

      if(!( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;
      if(!( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;


      ///// Acceptance Cut /////
      bool muplAcc = (
          ( (TMath::Abs(mupl_Reco->Eta()) <= 1.2) && (mupl_Reco->Pt() >=3.5) ) ||
          ( (TMath::Abs(mupl_Reco->Eta()) > 1.2)  && (TMath::Abs(mupl_Reco->Eta()) <= 2.1) && (mupl_Reco->Pt() >= 5.47-1.89*(TMath::Abs(mupl_Reco->Eta()))) ) ||
          ( (TMath::Abs(mupl_Reco->Eta()) > 2.1)  && (TMath::Abs(mupl_Reco->Eta()) <= 2.4) && (mupl_Reco->Pt() >= 1.5) ) 
          ) ;
      bool mumiAcc = (
          ( (TMath::Abs(mumi_Reco->Eta()) <= 1.2) && (mumi_Reco->Pt() >=3.5) ) ||
          ( (TMath::Abs(mumi_Reco->Eta()) > 1.2)  && (TMath::Abs(mumi_Reco->Eta()) <= 2.1) && (mumi_Reco->Pt() >= 5.47-1.89*(TMath::Abs(mumi_Reco->Eta()))) ) ||
          ( (TMath::Abs(mumi_Reco->Eta()) > 2.1)  && (TMath::Abs(mumi_Reco->Eta()) <= 2.4) && (mumi_Reco->Pt() >= 1.5) ) 
          );

      if ( !(muplAcc && mumiAcc) )
        continue ;

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

      if ( !(muplSoft && mumiSoft) ) 
        continue;   

      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) 
        continue;

      recoQQsign[irqq] = Reco_QQ_sign[irqq];     

      count++;     
      if(isMC){
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
          //        cout << "TnP ERROR !!!! ::: No matched L2 filter1 " << endl;
          continue;
        }
        bool mupl_L2Filter = ( (Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ? true : false ;
        bool mupl_L3Filter = ( (Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)) ) ? true : false ;
        bool mumi_L2Filter = ( (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ? true : false ;
        bool mumi_L3Filter = ( (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)) ) ? true : false ;
        //if(mupl_L2Filter == false || mumi_L2Filter == false){ continue;} 
        if(mupl_L2Filter == false || mumi_L2Filter == false){ cout << "TnP ERROR !!!! ::: No matched L2 filter2 " << endl; cout << endl; } 

        bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter) ? true : false;
        bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter) ? true : false;
        bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter) ? true : false;
        bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter) ? true : false;
        bool SelDone = false;
        /*
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
        //       if(SelDone == false || (tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1)){cout << "ERROR :: No muon filter combination selected !!!!" << endl; continue;}
        //	   { continue; }*/
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

        tnp_weight = tnp_weight * tnp_trig_weight_mupl * tnp_trig_weight_mumi;
        tnp_weight =1;
        //cout<<tnp_weight<<endl;
        counttnp++;
      }

      // Fill the output tree
      if ( JP_Reco->Eta() < 0 )  {  
        qxa[nDimu] = qx[HFp2];
        qya[nDimu] = qy[HFp2];
        qxb[nDimu] = qx[HFm2];
        qyb[nDimu] = qy[HFm2];

      }
      else {
        qxa[nDimu] = qx[HFm2];
        qya[nDimu] = qy[HFm2];
        qxb[nDimu] = qx[HFp2];
        qyb[nDimu] = qy[HFp2];
      }

      qxc[nDimu] = qx[trackmid2];
      qyc[nDimu] = qy[trackmid2];

      if(isMC) TnPweight[nDimu] = tnp_weight;
      mass[nDimu] = JP_Reco->M();
      pt[nDimu] = JP_Reco->Pt();
      pt1[nDimu] = mupl_Reco->Pt();
      pt2[nDimu] = mumi_Reco->Pt();
      y[nDimu] = JP_Reco->Rapidity();
      eta[nDimu] = JP_Reco->Eta();
      eta1[nDimu] = mupl_Reco->Eta();
      eta2[nDimu] = mumi_Reco->Eta();
      phi[nDimu] = JP_Reco->Phi();
      phi1[nDimu] = mupl_Reco->Phi();
      phi2[nDimu] = mumi_Reco->Phi();
      qxdimu[nDimu] = TMath::Cos(2*phi[nDimu]);
      qydimu[nDimu] = TMath::Sin(2*phi[nDimu]);
      qxmupl[nDimu] = TMath::Cos(2*phi1[nDimu]);
      qxmumi[nDimu] = TMath::Cos(2*phi2[nDimu]);
      qymupl[nDimu] = TMath::Sin(2*phi1[nDimu]);
      qymumi[nDimu] = TMath::Sin(2*phi2[nDimu]);
      ctau3D[nDimu] = Reco_QQ_ctau3D[irqq];
      ctau3DErr[nDimu] = Reco_QQ_ctauErr3D[irqq];
      ctau3Dtrue[nDimu] = Gen_QQ_ctau3D[Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]]];
      //cout<<ctau3D[nDimu]<<"-"<<ctau3Dtrue[nDimu]<<endl;
      //Ctau3DErrVar->setVal(ctau3DErr[nDimu]);
      //dataSet->add( *argSet);
      nDimu++;

    } // end of dimuon loop

    if(nDimu>0) mmevttree->Fill();  

  } //end of event loop

  cout<<"End fill Reco tree"<<endl;

  //////////////////////////////////////////////////////
  //////////////////////// Gen /////////////////////////
  //////////////////////////////////////////////////////

  cout<<""<<endl;
  cout<<"Start fill Gen tree"<<endl;

  for(int iev=0; iev<nevt ; ++iev)
  {
    //if(iev==1000)break;
    //if(iev==10000)break;
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
    nDimuGen = 0;

    for (Int_t irgen=0; irgen<Gen_QQ_size; ++irgen)
    {
      JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(irgen);
      mupl_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[irgen]);
      mumi_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[irgen]);

      evt = eventNb;

      count++;

      Genmass[nDimuGen] = JP_Gen->M();
      Genpt[nDimuGen] = JP_Gen->Pt();
      Geny[nDimuGen] = JP_Gen->Rapidity();
      Geneta[nDimuGen] = JP_Gen->Eta();
      Genphi[nDimuGen] = JP_Gen->Phi();
      Genpt1[nDimuGen] = mupl_Gen->Pt();
      Genpt2[nDimuGen] = mumi_Gen->Pt();
      Geneta1[nDimuGen] = mupl_Gen->Eta();
      Geneta2[nDimuGen] = mumi_Gen->Eta();
      Genphi1[nDimuGen] = mupl_Gen->Phi();
      Genphi2[nDimuGen] = mumi_Gen->Phi();
      Genctau3D[nDimuGen] = Gen_QQ_ctau3D[irgen];
      nDimuGen++;
    }

    if (nDimuGen>0 ) mmgentree->Fill();
  }

  //  mmtree->Write();  // Don't need to call Write() for trees
  cout << "count " << count << endl;
  cout << "counttnp " << counttnp << endl;
  newfile->cd();
  mmevttree->Write();
  mmgentree->Write();
  //datasSet->Write();
  newfile->Close();
} 
