#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBin.h"

using namespace std;

valErr getYield(float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0, int cLow=0, int cHigh=0);
void dndpt(int PR=1, int WRITE=1) {

  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1);
  gStyle->SetPadTickY(1);

  TH1::SetDefaultSumw2();

  //// modify by hand according to the pt range of the sample
  const int nPtBins=9;
  double ptBin[nPtBins+1]={3.0,4.0,5.0,6.0,7.0,8.0,10,12.0, 20.0, 50.0};
  const int nPtBinsMC=5;
  double ptBinMC[nPtBinsMC+1]={3.0,4.5,6.5,8.0,12.0,50.0};
  const int nYBins=6;
  double yBin[nYBins+1]={0.0,0.4,0.8,1.2,1.6,2.0,2.4};
  double frac[nPtBins]={0.1616, 0.18, 0.2, 0.21, 0.23, 0.24, 0.2985, 0.3740, 0.4138};
  // Get MC :
  float massLow = 2.6; float massHigh = 3.5;
  double ptMin = ptBinMC[0]; double ptMax = ptBinMC[nPtBinsMC];
  double yMin = yBin[0];     double yMax = yBin[nYBins];

  TH1D* hptData=new TH1D("hptData",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hptData1=new TH1D("hptData1",";p_{T}(GeV/c);",nPtBins,ptBin);
  //TH1D* hptData2=new TH1D("hptData2",";p_{T}(GeV/c);",nPtBins,ptBin);
  //TH1D* hptData3=new TH1D("hptData3",";p_{T}(GeV/c);",nPtBins,ptBin);
  //TH1D* hptData4=new TH1D("hptData4",";p_{T}(GeV/c);",nPtBins,ptBin);
  //TH1D* hptData5=new TH1D("hptData5",";p_{T}(GeV/c);",nPtBins,ptBin);
  //TH1D* hptData6=new TH1D("hptData6",";p_{T}(GeV/c);",nPtBins,ptBin);
  //TH1D* hptData7=new TH1D("hptData7",";|y|;",nYBins,yBin);

  TH1D *DataFit1 = new TH1D("DataFit1","; p_{T} (GeV/c) ; ", nPtBins, ptBin);

  TH1D *WeightFactor = new TH1D("WeightFactor","; p_{T} (GeV/c) ; ", nPtBinsMC, ptBinMC);

  TH1D *hptMC  = new TH1D("hptMC","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *hptMC1 = new TH1D("hptMC1","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *hptMC2 = new TH1D("hptMC2","; p_{T} (GeV/c) ; ", 100, 0, 50 );
  TH1D *hptMC3 = new TH1D("hptMC3","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *hptMC4 = new TH1D("hptMC4","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *hptMC5 = new TH1D("hptMC5","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *hptMC6 = new TH1D("hptMC6","; p_{T} (GeV/c) ; ",nPtBins,ptBin);
  TH1D *hptMC7 = new TH1D("hptMC7","; |y| ; ", nYBins, yBin);
  TH1D *MCfit  = new TH1D("MCfit","; p_{T} (GeV/c) ; ", nPtBinsMC, ptBinMC);

  TChain *tree = new TChain("mmepevt");
  TString f1;
  if(PR==0)  f1 ="../skimmedFiles/OniaFlowSkim_JpsiTrig_JPsi_isMC1_HFNom_201005.root";
  else if(PR==1) f1 ="../skimmedFiles/OniaFlowSkim_JpsiTrig_BtoJPsi_isMC1_HFNom_201005.root";
  tree->Add(f1.Data());

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
  float P[nMaxDimu];
  float Px[nMaxDimu];
  float Py[nMaxDimu];
  float Pz[nMaxDimu];
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

  TBranch *b_event;
  TBranch *b_cBin;
  TBranch *b_nDimu;
  TBranch *b_vz;
  TBranch *b_mass;
  TBranch *b_recoQQsign;
  TBranch *b_P;
  TBranch *b_Px;
  TBranch *b_Py;
  TBranch *b_Pz;
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
  tree -> SetBranchAddress("qxa", qxa, &b_qxa);
  tree -> SetBranchAddress("qxb", qxb, &b_qxb);
  tree -> SetBranchAddress("qxc", qxc, &b_qxc);
  tree -> SetBranchAddress("qxdimu", qxdimu, &b_qxdimu);
  tree -> SetBranchAddress("qya", qya, &b_qya);
  tree -> SetBranchAddress("qyb", qyb, &b_qyb);
  tree -> SetBranchAddress("qyc", qyc, &b_qyc);
  tree -> SetBranchAddress("qydimu", qydimu, &b_qydimu);

  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;
  for(int i=0; i<nEvt; i++){
    tree->GetEntry(i);

    for(int j=0; j<nDimu; j++){
      if ( !((fabs(eta1[j])>1.6) && (fabs(eta2[j])>1.6)))
        continue;
      if (  !( (mass[j] > massLow)
            && (mass[j] < massHigh)
            && ( pt[j] > ptMin)
            && ( pt[j] < ptMax)
            && ( fabs(y[j]) > 1.6)
            && ( fabs(y[j]) < 2.4) 
            && cBin>20
            && cBin<120 )
         )
        continue;
      hptMC->Fill      ( pt[j] );
      hptMC1->Fill     ( pt[j] );
      hptMC2->Fill     ( pt[j] );
    }
    //hptMC3->Fill     ( pt[i] );
    //hptMC6->Fill     ( pt[i] );
  }

  //------------------------------------------ Get Data :
  //FROM FINAL RESULTS
  for(int ipt=1;ipt<=nPtBins;ipt++)
  {
    valErr yieldAA;
    yieldAA = getYield(ptBin[ipt-1],ptBin[ipt],1.6,2.4,20,120);
    //yieldAA = getYield(8.0,12.0,0,2.4,20,120);
    if(PR==0){
    hptData->SetBinContent(ipt,yieldAA.val*(1-frac[ipt-1]));
    hptData->SetBinError(ipt,yieldAA.err);
    hptData1->SetBinContent(ipt,yieldAA.val*(1-frac[ipt-1]));
    hptData1->SetBinError(ipt,yieldAA.err);
    }
    if(PR==1){
    hptData->SetBinContent(ipt,yieldAA.val*frac[ipt-1]);
    hptData->SetBinError(ipt,yieldAA.err);
    hptData1->SetBinContent(ipt,yieldAA.val*frac[ipt-1]);
    hptData1->SetBinError(ipt,yieldAA.err);
    }
    //hptData2->SetBinContent(ipt,(yieldAA.val)+(yieldAA.err));
    //hptData3->SetBinContent(ipt,(yieldAA.val)-(yieldAA.err));
    //hptData4->SetBinContent(ipt,yieldAA.val);
    //hptData4->SetBinError(ipt,yieldAA.err);
    //hptData5->SetBinContent(ipt,yieldAA.val);
    //hptData5->SetBinError(ipt,yieldAA.err);
    //hptData6->SetBinContent(ipt,yieldAA.val);
    //hptData6->SetBinError(ipt,yieldAA.err);
    cout << Form("yield, pt  ") << ptBin[ipt-1] << " - " << ptBin[ipt] << " : " << yieldAA.val <<" +/- "<< yieldAA.err 
      << ", error : " << 100*yieldAA.err/yieldAA.val <<"%"<< ", Frac: "<<frac[ipt-1]<<endl;
  }
  ///////////////////////////////////Normalization///////////////////////////////////
  hptMC->Scale(1./hptMC->Integral());
  hptMC1->Scale(1./hptMC1->Integral());
  hptMC2->Scale(1./hptMC2->Integral());
  //hptMC3->Scale(1./hptMC3->Integral());
  //hptMC4->Scale(1./hptMC4->Integral());
  //hptMC5->Scale(1./hptMC5->Integral());
  //hptMC6->Scale(1./hptMC6->Integral());
  //hptMC7->Scale(1./hptMC7->Integral());

  hptData->Scale(1./hptData->Integral());
  hptData1->Scale(1./hptData1->Integral());
  //hptData2->Scale(1./hptData2->Integral());
  //hptData3->Scale(1./hptData3->Integral());
  //hptData4->Scale(1./hptData4->Integral());
  //hptData5->Scale(1./hptData5->Integral());
  //hptData6->Scale(1./hptData6->Integral());

  TH1ScaleByWidth(hptMC);
  TH1ScaleByWidth(hptData);

  handsomeTH1(hptMC,1);     
  handsomeTH1(hptData,1);   
  handsomeTH1(hptData1,1);   

  ///////////////////////////////////////////////////Fit Function///////////////////////////
  TF1* fitmc1;
  TF1* fitdata1;
  TF1* fitRatio1;

  fitmc1 = new TF1("fitmc1","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*x))/([0]*[1])),-[0]))",0,30);
  fitdata1 = new TF1("fitdata1","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*x))/([0]*[1])),-[0]))",0,30);
  //fitRatio1 = new TF1("fitRatio1","(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*x))/([0]*[1])),-[0]))/([2]-1)*([3]-2)/([2]*[3]*([2]*[3] + ([2]-2)*[2])) * x * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*x))/([2]*[3])),-[2])",0,30);
  fitRatio1 = new TF1("fitRatio1","( [0] + [1]*x + [2]*x*x +[4]*x*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",3.0,50);
  //fitRatio1 = new TF1("fitRatio1","TMath::Exp(-x/[0])*[1]",3.5,50);
  //fitRatio1->SetParameters(0.0005,4.5);
  //fitRatio1->SetParameters(0, 10);
  //fitRatio1->SetParameters(1, 10);

  TLegend *leg1 = new TLegend(0.65,0.75,0.85,0.85);
  leg1->AddEntry(hptData,"Data","p");
  leg1->AddEntry(hptMC,"MC","l");
  leg1->SetLineColor(kWhite);

  //TLegend *leg2 = new TLegend(0.5,0.65,0.87,0.85);
  //leg2->AddEntry(fitRatio1,Form("reweighting factor using fit(%s %dS)",fcollId.Data(), state),"l");
  //leg2->SetLineColor(kWhite);

  TCanvas* c_A =  new TCanvas("canvas_A","My plots",4,4,550,520);
  c_A->cd();
  TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0, 0.16, 0.98, 1.0);
  pad_A_1->SetTicks(1,1);
  pad_A_1->Draw(); pad_A_1->cd();
  //pad_A_1->Range(-2.857971,-0.0643152,55.88522,0.5337711);
  pad_A_1->SetFillColor(0);
  pad_A_1->SetBorderMode(0);
  pad_A_1->SetBorderSize(2);
  pad_A_1->SetTicks(1,1);
  pad_A_1->SetTopMargin(0.05646528);
  pad_A_1->SetFrameBorderMode(0);
  pad_A_1->SetFrameBorderMode(0);
  hptData->Draw();
  hptMC->Draw("same hist");
  leg1->Draw("same");
  hptData->SetAxisRange(0.,hptData->GetMaximum()+0.05,"Y");
  hptData->GetXaxis()->SetLabelSize(0);
  hptData->GetXaxis()->SetTitleSize(0);
  //TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2",0,0.09,0.98,0.23);
  TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0, 0.006, 0.98, 0.227);
  c_A->cd();
  pad_A_2->Draw();
  pad_A_2->cd();
  //pad_A_2->Range(-2.857971,-3.915054,55.88522,5.062366);
  pad_A_2->SetFillColor(0);
  pad_A_2->SetBorderMode(0);
  pad_A_2->SetBorderSize(2);
  pad_A_2->SetTicks(1,1);
  //pad_A_2->SetTopMargin(0.00694694);
  pad_A_2->SetBottomMargin(0.4361001);
  //pad_A_2->SetFrameBorderMode(0);
  //pad_A_2->SetFrameBorderMode(0);
  hptData1->Divide(hptMC1);
  //hptData2->Divide(hptMC2);
  //hptData3->Divide(hptMC3);
  //hptData1->Fit(fitRatio1,"IE","",0,30);
  //TFitResultPtr r = hptData1->Fit(fitRatio1,"S");
  //r.Get()->Print("V");
  //fitRatio2->FixParameter(0,fitRatio1->GetParameter(0)+fitRatio1->GetParError(0));
  //fitRatio2->FixParameter(1,fitRatio1->GetParameter(1)+fitRatio1->GetParError(1));
  //fitRatio2->FixParameter(2,fitRatio1->GetParameter(2)+fitRatio1->GetParError(2));
  //fitRatio2->FixParameter(3,fitRatio1->GetParameter(3)+fitRatio1->GetParError(3));
  //fitRatio3->FixParameter(0,fitRatio1->GetParameter(0)-fitRatio1->GetParError(0));
  //fitRatio3->FixParameter(1,fitRatio1->GetParameter(1)-fitRatio1->GetParError(1));
  //fitRatio3->FixParameter(2,fitRatio1->GetParameter(2)-fitRatio1->GetParError(2));
  //fitRatio3->FixParameter(3,fitRatio1->GetParameter(3)-fitRatio1->GetParError(3));
  //cout<<"Fit Parameter 1"<<endl;
  //cout<<fitRatio1->GetParameter(0)<<", "<<fitRatio1->GetParameter(1)<<", "<<fitRatio1->GetParameter(2)<<", "<<fitRatio1->GetParameter(3)<<endl;
  //cout<<"Fit Parameter 2"<<endl;
  //cout<<fitRatio2->GetParameter(0)<<", "<<fitRatio2->GetParameter(1)<<", "<<fitRatio2->GetParameter(2)<<", "<<fitRatio2->GetParameter(3)<<endl;
  //cout<<"Fit Parameter 3"<<endl;
  //cout<<fitRatio3->GetParameter(0)<<", "<<fitRatio3->GetParameter(1)<<", "<<fitRatio3->GetParameter(2)<<", "<<fitRatio3->GetParameter(3)<<endl;
  hptData1->Draw();
  hptData1->GetXaxis()->SetTitleOffset(1.2) ;
  hptData1->GetXaxis()->SetTitleSize(0.15) ;
  hptData1->GetXaxis()->CenterTitle();
  hptData1->GetXaxis()->SetLabelOffset(0.04) ;
  hptData1->GetXaxis()->SetLabelSize(0.15) ;
  hptData1->GetXaxis()->SetTickSize(0.03);
  hptData1->GetYaxis()->SetTickSize(0.04);
  hptData1->GetYaxis()->SetNdivisions(404);
  hptData1->GetYaxis()->SetTitle("Data/MC");
  hptData1->GetYaxis()->SetTitleOffset(0.25) ;
  hptData1->GetYaxis()->SetTitleSize(0.15) ;
  hptData1->GetYaxis()->CenterTitle();
  hptData1->GetYaxis()->SetLabelSize(0.15) ;
  hptData1->GetYaxis()->SetTickSize(0.04);
  hptData1->GetYaxis()->SetNdivisions(404);
  //hptData2->SetLineColor(kGreen+2);
  //hptData3->SetLineColor(kRed+2);
  //fitRatio1->SetLineColor(kGreen+2);
  //fitRatio2->SetLineColor(kRed+2);
  //fitRatio3->SetLineColor(kBlue+2);
  //fitRatio1->Draw("same");
  //fitRatio2->Draw("same");
  //fitRatio3->Draw("same");
  hptData1->SetAxisRange(-0.5,5.,"Y");
  hptData1->Fit(fitRatio1,"IE","",3.0,50);
  TFitResultPtr r = hptData1->Fit(fitRatio1,"S","",3.0,50);
  //TFitResultPtr r = hptData1->Fit("expo","S");
  r.Get()->Print("V");
  //leg2->Draw("same");
  jumSun(ptMin,1,ptMax,1);

  TCanvas* c2 =  new TCanvas("c2","",604, 0, 800, 600);
  hptMC2->Draw("same hist");
  ////fitRatio1->Draw();
  ////  TCanvas* c1 =  new TCanvas("c1","",1200, 800);
  ////  c1->Divide(3,2);
  ////  //MC plot
  ////  c1->cd(1) ;
  ////  //hptMC->Fit(fitmc1,"I");
  ////  hptMC->Draw();
  ////  gPad->SetLogy();
  ////  //Data plot
  ////  cout<<"FIND!!!!!"<<endl;
  ////  c1->cd(2);
  ////  gPad->SetLogy();
  ////  hptData->Fit(fitdata1,"I");
  ////  hptData->Draw();
  ////  fitdata2->FixParameter(0,fitdata1->GetParameter(0));
  ////  fitdata2->FixParameter(1,fitdata1->GetParameter(1));
  ////  fitdata2->FixParameter(2,fitdata1->GetParameter(2));
  ////  fitdata3->FixParameter(0,fitdata1->GetParameter(0));
  ////  fitdata3->FixParameter(1,fitdata1->GetParameter(1));
  ////  fitdata3->FixParameter(2,fitdata1->GetParameter(2));
  ////  fitdata2->Draw("same");
  ////  fitdata3->Draw("same");
  ////  //hptData->SetAxisRange(0.,0.3,"Y");
  ////  //Weighting Factor
  ////  //handsomeTH1(WeightFactor1,1);   handsomeTH1(WeightFactor2,1); handsomeTH1(WeightFactor3,1);
  ////
  ////  c1->cd(3);
  ////  cout<<"FIND!!!!!"<<endl;
  ////  float fittmp1, fittmp2, fittmp3;
  ////  float mctmp1, mctmp2, mctmp3;
  ////  for(int i =0; i<nPtBinsMC; ++i){
  ////    fittmp1=fitdata1->Eval(hptMC1->GetBinCenter(i+1));
  ////    fittmp2=fitdata2->Eval(hptMC1->GetBinCenter(i+1));
  ////    fittmp3=fitdata3->Eval(hptMC1->GetBinCenter(i+1));
  ////    //mctmp1=fitmc1->Eval(hptMC1->GetBinCenter(i+1));
  ////    //fittmp1=fitdata1->Integral(ptBinMC[i],ptBinMC[i+1]);
  ////    //fittmp2=fitdata2->Integral(ptBinMC[i],ptBinMC[i+1]);
  ////    //fittmp3=fitdata3->Integral(ptBinMC[i],ptBinMC[i+1]);
  ////    WeightFactor1->SetBinContent(i+1,fittmp1);
  ////    WeightFactor2->SetBinContent(i+1,fittmp2);
  ////    WeightFactor3->SetBinContent(i+1,fittmp3);
  ////    //WeightDeno1->SetBinContent(i+1,mctmp1);
  ////  }
  ////  //WeightFactor1->Divide(WeightDeno1);
  ////  WeightFactor1->Divide(hptMC1);
  ////  WeightFactor2->Divide(hptMC2);
  ////  WeightFactor3->Divide(hptMC3);
  ////  WeightFactor1->SetMarkerColor(kBlack);     WeightFactor2->SetMarkerColor(kRed+1);     WeightFactor3->SetMarkerColor(kGreen+2);
  ////  WeightFactor1->SetMarkerStyle(1);          WeightFactor2->SetMarkerStyle(2);          WeightFactor3->SetMarkerStyle(2);
  ////  WeightFactor1->Draw();
  ////  WeightFactor2->Draw("same");
  ////  WeightFactor3->Draw("same");
  ////  WeightFactor1->SetAxisRange(0.,2.,"Y");
  ////  if(state==2){
  ////    WeightFactor1->SetAxisRange(0.,6.,"Y");}
  ////  //hptData6->Divide(hptMC6);
  ////  //hptData6->Draw();
  ////  jumSun(0,1,30,1);
  ////  //for(int i =0; i<nPtBins; ++i){
  ////  //WeightFactor1->SetBinContent( i+1, hptData6->GetBinContent(i+1));
  ////  //}
  ////  //WeightFactor1->Draw("same");
  ////  //Fit-Data Ratio
  ////  c1->cd(4);
  ////  hptMC4->Draw("hist");
  ////  hptData4->Draw("same");
  ////  hptMC4->SetAxisRange(0.,0.165,"Y");
  ////  jumSun(0,1,30,1);
  ////
  ////  c1->cd(5);
  ////  hptMC5->Draw("hist");
  ////  //hptData4->Draw("same");
  ////  hptData5->Draw("same");
  ////  hptMC5->SetAxisRange(0.,0.165,"Y");
  ////  jumSun(0,1,30,1);
  ////
  ////  c1->cd(6);
  ////  float fitval1, fitval2, fitval3;
  ////  float fitint1, fitint2, fitint3;
  ////  handsomeTH1(DataFit1,1);  handsomeTH1(DataFit2,1); handsomeTH1(DataFit3,1);
  ////  for(int i =0; i<nPtBins; i++){
  ////    fitval1=fitdata1->Eval(hptData1->GetBinCenter(i+1));
  ////    fitval2=fitdata2->Eval(hptData2->GetBinCenter(i+1));
  ////    fitval3=fitdata3->Eval(hptData3->GetBinCenter(i+1));
  ////    fitint1=fitdata1->Integral(ptBin[i],ptBin[i+1]);
  ////    fitint2=fitdata2->Integral(ptBin[i],ptBin[i+1]);
  ////    fitint3=fitdata3->Integral(ptBin[i],ptBin[i+1]);
  ////    DataFit1->SetBinContent(i+1,fitint1);
  ////    DataFit2->SetBinContent(i+1,fitint2);
  ////    DataFit3->SetBinContent(i+1,fitint3);
  ////  }
  ////  TH1ScaleByWidth(DataFit1);  TH1ScaleByWidth(DataFit2);  TH1ScaleByWidth(DataFit3);
  ////  for(int i=1; i<=nPtBins; i++){
  ////    cout<<"  ["<<ptBin[i-1]<<"-"<<ptBin[i]<<"] Data :"<<hptData1->GetBinContent(i)<<" Fit Int : "<<DataFit1->GetBinContent(i)
  ////      <<" , Ratio : "<<DataFit1->GetBinContent(i)/hptData1->GetBinContent(i)<<endl;
  ////    //", Fit val : "<<fitval1<<endl;
  ////  }
  ////  DataFit1->SetMarkerColor(kBlack); DataFit2->SetMarkerColor(kRed+1); DataFit3->SetMarkerColor(kGreen+2);
  ////  DataFit1->SetLineColor(kBlack); DataFit2->SetLineColor(kRed+1); DataFit3->SetLineColor(kGreen+2);
  ////  DataFit1->Divide(hptData1);
  ////  DataFit2->Divide(hptData2);
  ////  DataFit3->Divide(hptData3);
  ////  DataFit1->Draw();
  ////  DataFit2->Draw("same");
  ////  DataFit3->Draw("same");
  ////  DataFit1->SetAxisRange(0.,2.,"Y");
  ////  jumSun(0,1,30,1);
  ////
  //c1->SaveAs(Form("dNdpt_plot_%s_%dS.pdf",fcollId.Data(),state));
  //c1->SaveAs(Form("dNdpt_plot_%s_%dS.png",fcollId.Data(),state));

  if(WRITE==1){
    TFile *fJpsipb = new TFile("WeightedFcN_fit/ratioDataMC_AA_btojpsi_DATA_1s.root","RECREATE");
    fJpsipb->cd();
    hptData1->SetName("WeightFactor");
    hptData1->Write();
    fitRatio1->SetName("dataMC_Ratio1");
    fitRatio1->Write();
  }
}

//Get Yield
valErr getYield(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh) {
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);
  TFile* inf = new TFile(Form("../roots/2DFit_1207/MassFitResult_Inclusive_%s.root", kineLabel.Data()));
  //RooWorkspace* ws = (RooWorkspace*)inf->Get("workspace");
  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  valErr ret;
  ret.val = fitResults->GetBinContent(1);
  ret.err = fitResults->GetBinError(1);
  //cout << kineLabel << ": " << " & " << ret.val << " $\pm$ " << ret.err << " & " <<ws->var("nBkg")->getVal() << " $\pm$ "<< ws->var("nBkg")->getError() << "\\\\" << endl;
  return ret;
}

