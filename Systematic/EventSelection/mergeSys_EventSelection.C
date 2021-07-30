#include "../../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../../cutsAndBin.h"
#include "../../Style.h"
using namespace std;

TLegend *leg = new TLegend(0.55,0.2, 0.85,0.4,NULL,"brNDC");

void mergeTwoInMax( TH1D* h0, TH1D* h1, TH1D* h2, const char* str);

void mergeSys_EventSelection() { 
  
  TH1::SetDefaultSumw2();

  TH1D*h_Pt1_PR[10];
  TH1D*h_Pt1_NP[10];
 
  // 1 : Event Selection Up
  TFile* f1 = new TFile("EventSelection_HFUp_sys.root");
  h_Pt1_PR[1] = (TH1D*) f1->Get("hPtSys_PR");
  h_Pt1_NP[1] = (TH1D*) f1->Get("hPtSys_NP");
  
  // 2 : Event Selection Down
  TFile* f2 = new TFile("EventSelection_HFDown_sys.root");
  h_Pt1_PR[2] = (TH1D*) f2->Get("hPtSys_PR");
  h_Pt1_NP[2] = (TH1D*) f2->Get("hPtSys_NP");

  //const char* str_ptPR[7] = {
  //  "p_{T} 3-4.5 GeV", "p_{T} 4.5-6.5 GeV", "p_{T} 6.5-7.5 GeV", "p_{T} 7.5-9 GeV", "p_{T} 9-12 GeV", "p_{T} 12-15 GeV", "p_{T} 12-50 GeV"};
  //const char* str_ptNP[5] = {
  //  "p_{T} 3-6.5 GeV", "p_{T} 6.5-9 GeV", "p_{T} 9-12 GeV", "p_{T} 12-15 GeV", "p_{T} 12-50 GeV"};
  const char* str_ptPR[7] = {
	"p_{T} 3-6.5 GeV/c, 1.6 < |y| < 2.4 p_{T} 6.5-50 GeV/c, |y| < 2.4"};
  const char* str_ptNP[5] = {
	"p_{T} 3-6.5 GeV", "p_{T} 6.5-9 GeV", "p_{T} 9-12 GeV", "p_{T} 12-15 GeV", "p_{T} 12-50 GeV"};

  h_Pt1_PR[0]  = (TH1D*) h_Pt1_PR[1]->Clone("hPtSys_PR"); h_Pt1_PR[0]->Reset();
  h_Pt1_NP[0]  = (TH1D*) h_Pt1_NP[1]->Clone("hPtSys_NP"); h_Pt1_NP[0]->Reset();

  mergeTwoInMax(h_Pt1_PR[0],h_Pt1_PR[1],h_Pt1_PR[2],str_ptPR[0]);
  mergeTwoInMax(h_Pt1_NP[0],h_Pt1_NP[1],h_Pt1_NP[2],str_ptNP[0]);

  TFile* fout = new TFile("EventSelection_sys.root","recreate");
  fout->cd();
  h_Pt1_PR[0]->Write();
  h_Pt1_NP[0]->Write();

  fout->Close();

}

void mergeTwoInMax( TH1D* h0, TH1D* h1, TH1D* h2, const char* str) {
  gStyle->SetOptStat(0);
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a0 = max(a1,a2);
    cout<<a0<<endl;
    h0->SetBinContent( i, a0);
  } 

  //TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);

  //h0->SetAxisRange(0,0.1,"Y");
  //h0->SetYTitle("Difference Uncertainty");
  //handsomeTH1(h0,        1); h0->SetLineWidth(2); h0->DrawCopy("hist");
  //handsomeTH1(h1,        2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  //handsomeTH1(h2,        3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  //
  //TLegend *leg1 = new TLegend(0.55,0.6, 0.85,0.9,NULL,"brNDC");
  //easyLeg(leg1,Form("%s, J#varPsi",h0->GetName()));
  //leg1->AddEntry(h0,"Total","l");
  //leg1->AddEntry(h1,"Signal PDF","l");
  //leg1->AddEntry(h2,"Background PDF","l");
  //leg1->Draw();
  //
  //TLatex* globtex = new TLatex();
  //globtex->SetNDC();
  //globtex->SetTextAlign(12); //left-center
  //globtex->SetTextFont(42);
  //globtex->SetTextSize(0.035);
  //globtex->DrawLatex(0.15, 0.87, Form("%s",str));
  //c0->SaveAs(Form("pdfFiles/Unc_%s.pdf", h0->GetName() ));
}
