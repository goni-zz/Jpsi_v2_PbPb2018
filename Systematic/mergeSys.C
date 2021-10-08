#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBin.h"
#include "../Style.h"
using namespace std;

TLegend *leg = new TLegend(0.55,0.2, 0.85,0.4,NULL,"brNDC");
void mergeEightInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D*h5=0, TH1D*h6=0, TH1D*h7=0, TH1D*h8=0, const char* str = " ");
void mergeTwoInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, const char* str = " ");
void mergeThreeInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, const char* str = " ");
void mergeFourInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, const char* str = " ");
void mergeFiveInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D* h5=0, const char* str = " ");
void mergeSixInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D* h5=0, TH1D* h6=0, TH1D* h7=0, TH1D* h8=0, TH1D* h9=0, const char* str = " ");

void mergeSys(int nbin=1) { 
  
  TH1::SetDefaultSumw2();

  TH1D*h_Pt1_PR[10];
  TH1D*h_Pt1_NP[10];
 
  TString bin;
  if(nbin==1) bin="10_60";
  else if(nbin==2) bin="20_40";
  else if(nbin==3) bin="0_180";
  // 1 : signal PDF var
  TFile* f1 = new TFile("SignalPDFvariation/SignalPDFVariation_sys.root");
  h_Pt1_PR[1] = (TH1D*) f1->Get("hPtSys_PR");
  h_Pt1_NP[1] = (TH1D*) f1->Get("hPtSys_NP");


  // 2 : background PDF var
  TFile* f2 = new TFile(Form("BkgPDFvariation/BkgPDFVariation_sys_%s.root",bin.Data()));
  h_Pt1_PR[2] = (TH1D*) f2->Get("hPtSys_PR");
  h_Pt1_NP[2] = (TH1D*) f2->Get("hPtSys_NP");

  // 3 : efficiency
  TFile* f3 = new TFile(Form("Efficiency_PtW/Efficiency_PtW_sys_%s.root",bin.Data()));
  h_Pt1_PR[3] = (TH1D*) f3->Get("hPtSys_PR");
  h_Pt1_NP[3] = (TH1D*) f3->Get("hPtSys_NP");

  // 4 : acceptance
  TFile* f4 = new TFile(Form("Acceptance/Acceptance_sys_%s.root",bin.Data()));
  h_Pt1_PR[4] = (TH1D*) f4->Get("hPtSys_PR");
  h_Pt1_NP[4] = (TH1D*) f4->Get("hPtSys_NP");

  // 5 : TnP Correction
  TFile* f5 = new TFile("Efficiency_TnP/Efficiency_TnP_sys.root");
  h_Pt1_PR[5] = (TH1D*) f5->Get("hPtSys_PR");
  h_Pt1_NP[5] = (TH1D*) f5->Get("hPtSys_NP");

  // 6 : v2 background func. var
  TFile* f6 = new TFile(Form("v2Bkg/v2Bkg_sys_%s.root",bin.Data()));
  h_Pt1_PR[6] = (TH1D*) f6->Get("hPtSys_PR");
  h_Pt1_NP[6] = (TH1D*) f6->Get("hPtSys_NP");

  // 7 : Event Selection
  TFile* f7 = new TFile(Form("EventSelection/EventSelection_sys_%s.root",bin.Data()));
  h_Pt1_PR[7] = (TH1D*) f7->Get("hPtSys_PR");
  h_Pt1_NP[7] = (TH1D*) f7->Get("hPtSys_NP");

  // 8 : b-fraction 
  TFile* f8 = new TFile("b_fraction/bfraction_sys.root");
  h_Pt1_PR[8] = (TH1D*) f8->Get("hPtSys_PR");
  h_Pt1_NP[8] = (TH1D*) f8->Get("hPtSys_NP");

  // 9 : SignalParVaraition 
  TFile* f9 = new TFile(Form("SignalParams/SignalParVariation_sys_%s.root",bin.Data()));
  h_Pt1_PR[9] = (TH1D*) f9->Get("hPtSys_PR");
  h_Pt1_NP[9] = (TH1D*) f9->Get("hPtSys_NP");
  
  //const char* str_ptPR[7] = {
  //  "p_{T} 3-4.5 GeV", "p_{T} 4.5-6.5 GeV", "p_{T} 6.5-7.5 GeV", "p_{T} 7.5-9 GeV", "p_{T} 9-12 GeV", "p_{T} 12-15 GeV", "p_{T} 12-50 GeV"};
  //const char* str_ptNP[5] = {
  //  "p_{T} 3-6.5 GeV", "p_{T} 6.5-9 GeV", "p_{T} 9-12 GeV", "p_{T} 12-15 GeV", "p_{T} 12-50 GeV"};
  const char* str_ptPR[7] = {
	Form("Centrality_%s%",bin.Data())};
  const char* str_ptNP[5] = {
	Form("Centrality_%s%",bin.Data())};

  h_Pt1_PR[0]  = (TH1D*) h_Pt1_PR[1]->Clone("hpt_PR_merged"); h_Pt1_PR[0]->Reset();
  h_Pt1_NP[0]  = (TH1D*) h_Pt1_NP[1]->Clone("hpt_NP_merged"); h_Pt1_NP[0]->Reset();

  //mergeEightInQuad(str_int);  
  //mergeTwoInQuad(h_Pt1_PR[0],h_Pt1_PR[1],h_Pt1_PR[2],str_ptPR[0]);
  //mergeThreeInQuad(h_Pt1_PR[0],h_Pt1_PR[1],h_Pt1_PR[2],h_Pt1_PR[4],str_ptPR[0]);
  //mergeFourInQuad(h_Pt1_PR[0],h_Pt1_PR[1],h_Pt1_PR[2],h_Pt1_PR[3],h_Pt1_PR[4],str_ptPR[0]);
  //mergeFiveInQuad(h_Pt1_PR[0],h_Pt1_PR[1],h_Pt1_PR[2],h_Pt1_PR[3],h_Pt1_PR[4],h_Pt1_NP[7],str_ptPR[0]);
  //mergeFiveInQuad(h_Pt1_NP[0],h_Pt1_NP[1],h_Pt1_NP[2],h_Pt1_NP[3],h_Pt1_NP[4],h_Pt1_NP[7],str_ptNP[0]);
  //mergeSixInQuad(h_Pt1_PR[0],h_Pt1_PR[1],h_Pt1_PR[2],h_Pt1_PR[3],h_Pt1_PR[4],h_Pt1_PR[5],h_Pt1_PR[6],h_Pt1_PR[7],h_Pt1_PR[8],str_ptPR[0]);
  //mergeSixInQuad(h_Pt1_NP[0],h_Pt1_NP[1],h_Pt1_NP[2],h_Pt1_NP[3],h_Pt1_NP[4],h_Pt1_NP[5],h_Pt1_NP[6],h_Pt1_NP[7],h_Pt1_NP[8],str_ptNP[0]);
  mergeSixInQuad(h_Pt1_PR[0],h_Pt1_PR[1],h_Pt1_PR[2],h_Pt1_PR[3],h_Pt1_PR[4],h_Pt1_PR[5],h_Pt1_PR[6],h_Pt1_PR[7],h_Pt1_PR[8],h_Pt1_PR[9],str_ptPR[0]);
  mergeSixInQuad(h_Pt1_NP[0],h_Pt1_NP[1],h_Pt1_NP[2],h_Pt1_NP[3],h_Pt1_NP[4],h_Pt1_NP[5],h_Pt1_NP[6],h_Pt1_NP[7],h_Pt1_NP[8],h_Pt1_NP[9],str_ptNP[0]);
    
  TFile* fout = new TFile(Form("merged_sys_%s.root",bin.Data()),"recreate");
  fout->cd();
  h_Pt1_PR[0]->Write();
  h_Pt1_NP[0]->Write();

  fout->Close();

}

void mergeEightInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7, TH1D* h8, const char* str) {
  gStyle->SetOptStat(0);
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a5 = h5->GetBinContent(i);
    float a6 = h6->GetBinContent(i);
    float a7 = h7->GetBinContent(i);
    float a8 = h8->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 + a5*a5 + a6*a6 + a7*a7 + a8*a8);
    h0->SetBinContent( i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);

  h0->SetAxisRange(0,0.1,"Y");
  h0->SetYTitle("Difference Uncertainty");
  handsomeTH1(h0,        1); h0->SetLineWidth(2); h0->DrawCopy("hist");
  handsomeTH1(h1,        2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2,        3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(h3,        4); h3->SetLineWidth(2); h3->DrawCopy("hist same");
  handsomeTH1(h4,        5); h4->SetLineWidth(2); h4->DrawCopy("hist same");
  handsomeTH1(h5,        6); h5->SetLineWidth(2); h5->DrawCopy("hist same");
  handsomeTH1(h6,        7); h6->SetLineWidth(2); h6->DrawCopy("hist same");
  handsomeTH1(h7,        14); h7->SetLineWidth(2); h7->DrawCopy("hist same");
  handsomeTH1(h8,        20); h8->SetLineWidth(2); h8->DrawCopy("hist same");
  
  TLegend *leg1 = new TLegend(0.55,0.6, 0.85,0.9,NULL,"brNDC");
  easyLeg(leg1,"Merged Uncertainty");
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"Signal PDF","l");
  leg1->AddEntry(h2,"Background PDF","l");
  leg1->AddEntry(h3,"efficiency","l");
  leg1->AddEntry(h4,"Acceptance","l");
  leg1->AddEntry(h5,"b-fraction","l");
  leg1->AddEntry(h6,"v_{2} background PDF","l");
  leg1->AddEntry(h7,"Event Selection","l");
  leg1->AddEntry(h8,"Tag-And-Probe","l");
  leg1->Draw();
  
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.035);
  globtex->DrawLatex(0.15, 0.87, Form("%s",str));
  c0->SaveAs(Form("pdfFiles/Unc_%s.pdf", h0->GetName() ));
}

void mergeTwoInQuad( TH1D* h0, TH1D* h1, TH1D* h2, const char* str) {
  gStyle->SetOptStat(0);
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2);
    h0->SetBinContent( i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);

  h0->SetAxisRange(0,0.1,"Y");
  h0->SetYTitle("Difference Uncertainty");
  handsomeTH1(h0,        1); h0->SetLineWidth(2); h0->DrawCopy("hist");
  handsomeTH1(h1,        2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2,        3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  
  TLegend *leg1 = new TLegend(0.55,0.6, 0.9,0.9,NULL,"brNDC");
  //easyLeg(leg1,Form("%s, J/#psi",h0->GetName()));
  easyLeg(leg1,"Merged Uncertainty");
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"Signal PDF","l");
  leg1->AddEntry(h2,"Background PDF","l");
  leg1->Draw();
  
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.035);
  globtex->DrawLatex(0.15, 0.87, Form("%s",str));
  c0->SaveAs(Form("pdfFiles/Unc_%s.pdf", h0->GetName() ));
}

void mergeThreeInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, const char* str) {
  gStyle->SetOptStat(0);
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3);
    h0->SetBinContent(i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);

  h0->SetAxisRange(0,0.1,"Y");
  h0->SetYTitle("Difference Uncertainty");
  handsomeTH1(h0, 1); h0->SetLineWidth(2); h0->DrawCopy("hist");
  handsomeTH1(h1, 2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2, 3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(h3, 4); h3->SetLineWidth(2); h3->DrawCopy("hist same");
  
  TLegend *leg1 = new TLegend(0.55,0.6, 0.9,0.9,NULL,"brNDC");
  //easyLeg(leg1,Form("%s, J/#psi",h0->GetName()));
  easyLeg(leg1,"Merged Uncertainty");
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"Signal PDF","l");
  leg1->AddEntry(h2,"Background PDF","l");
  leg1->AddEntry(h3,"Acceptance","l");
  leg1->Draw();
  
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.035);
  globtex->DrawLatex(0.15, 0.87, Form("%s",str));
  c0->SaveAs(Form("pdfFiles/Unc_%s.pdf", h0->GetName() ));
}

void mergeFourInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, const char* str) {
  gStyle->SetOptStat(0);
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 +a4*a4);
    h0->SetBinContent(i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);

  h0->SetAxisRange(0,0.1,"Y");
  h0->SetYTitle("Difference Uncertainty");
  handsomeTH1(h0, 1); h0->SetLineWidth(2); h0->DrawCopy("hist");
  handsomeTH1(h1, 2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2, 3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(h3, 4); h3->SetLineWidth(2); h3->DrawCopy("hist same");
  handsomeTH1(h4, 5); h4->SetLineWidth(2); h4->DrawCopy("hist same");
  
  TLegend *leg1 = new TLegend(0.55,0.6, 0.9,0.9,NULL,"brNDC");
  //easyLeg(leg1,Form("%s, J/#psi",h0->GetName()));
  easyLeg(leg1,"Merged Uncertainty");
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"Signal PDF","l");
  leg1->AddEntry(h2,"Background PDF","l");
  leg1->AddEntry(h3,"Efficiency","l");
  leg1->AddEntry(h4,"Acceptance","l");
  leg1->Draw();
  
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.035);
  globtex->DrawLatex(0.15, 0.87, Form("%s",str));
  c0->SaveAs(Form("pdfFiles/Unc_%s.pdf", h0->GetName() ));
}

void mergeFiveInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, const char* str) {
  gStyle->SetOptStat(0);
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a5 = h5->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 + a5*a5);
    h0->SetBinContent(i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);

  h0->SetAxisRange(0,0.1,"Y");
  h0->SetYTitle("Difference Uncertainty");
  handsomeTH1(h0, 1); h0->SetLineWidth(2); h0->DrawCopy("hist");
  handsomeTH1(h1, 2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2, 3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(h3, 4); h3->SetLineWidth(2); h3->DrawCopy("hist same");
  handsomeTH1(h4, 5); h4->SetLineWidth(2); h4->DrawCopy("hist same");
  handsomeTH1(h5, 6); h5->SetLineWidth(2); h5->DrawCopy("hist same");
  
  TLegend *leg1 = new TLegend(0.55,0.6, 0.9,0.9,NULL,"brNDC");
  //easyLeg(leg1,Form("%s, J/#psi",h0->GetName()));
  easyLeg(leg1,"Merged Uncertainty");
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"Signal PDF","l");
  leg1->AddEntry(h2,"Background PDF","l");
  leg1->AddEntry(h3,"Efficiency","l");
  leg1->AddEntry(h4,"Acceptance","l");
  leg1->AddEntry(h5,"Event Selection","l");
  leg1->Draw();
  
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.035);
  globtex->DrawLatex(0.15, 0.87, Form("%s",str));
  c0->SaveAs(Form("pdfFiles/Unc_%s.pdf", h0->GetName() ));
}

void mergeSixInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7, TH1D* h8, TH1D* h9, const char* str) {
  gStyle->SetOptStat(0);
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a5 = h5->GetBinContent(i);
    float a6 = h6->GetBinContent(i);
    float a7 = h7->GetBinContent(i);
    float a8 = h8->GetBinContent(i);
    float a9 = h9->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 + a5*a5 + a6*a6 +a7*a7 +a8*a8 +a9*a9);
    h0->SetBinContent(i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);

  h0->SetAxisRange(0,0.04,"Y");
  h0->SetYTitle("Difference Uncertainty");
  handsomeTH1(h0, 1); h0->SetLineWidth(2); h0->DrawCopy("hist");
  handsomeTH1(h1, 2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2, 3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(h3, 4); h3->SetLineWidth(2); h3->DrawCopy("hist same");
  handsomeTH1(h4, 5); h4->SetLineWidth(2); h4->DrawCopy("hist same");
  handsomeTH1(h5, 6); h5->SetLineWidth(2); h5->DrawCopy("hist same");
  handsomeTH1(h6, 7); h6->SetLineWidth(2); h6->DrawCopy("hist same");
  handsomeTH1(h7, 11); h7->SetLineWidth(2); h7->DrawCopy("hist same");
  handsomeTH1(h8, 12); h8->SetLineWidth(2); h8->DrawCopy("hist same");
  handsomeTH1(h9, 13); h8->SetLineWidth(2); h9->DrawCopy("hist same");
  
  TLegend *leg1 = new TLegend(0.55,0.6, 0.9,0.9,NULL,"brNDC");
  //easyLeg(leg1,Form("%s, J/#psi",h0->GetName()));
  easyLeg(leg1,"Merged Uncertainty");
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"Signal PDF","l");
  leg1->AddEntry(h2,"Background PDF","l");
  leg1->AddEntry(h3,"Efficiency","l");
  leg1->AddEntry(h4,"Acceptance","l");
  leg1->AddEntry(h5,"Tag and Probe","l");
  leg1->AddEntry(h6,"v_{2} Background","l");
  leg1->AddEntry(h7,"Event Selection","l");
  leg1->AddEntry(h8,"b_fraction","l");
  leg1->AddEntry(h9,"Signal Par","l");
  leg1->Draw();
  
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.035);
  globtex->DrawLatex(0.15, 0.87, Form("%s",str));
  c0->SaveAs(Form("pdfFiles/Unc_%s_%s.pdf", h0->GetName(), str));
}

void mergeThreeInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3) {
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 );
    h0->SetBinContent( i, a0);
  }
}

void DrawSysMerged(TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3,TH1D* h4,TH1D* h5,TH1D* h6) {
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 );
    h0->SetBinContent( i, a0);
  }
}
