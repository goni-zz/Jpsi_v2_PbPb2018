#include "../../tdrstyle.C"
#include "../../CMS_lumi_v2mass.C"

void drawJpsi_npr_v2_pt_7bin(){
  //gROOT->Macro("~/rootlogon.C");
  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;

  gStyle->SetEndErrorSize(3);
  gStyle->SetPadTopMargin(0.07);

  TFile *fPt1 = new TFile("v2_D0method_pt3.0-4.5_y1.6-2.4_muPt0.0_centrality20-120.root","read");
  TFile *fPt2 = new TFile("v2_D0method_pt4.5-6.5_y1.6-2.4_muPt0.0_centrality20-120.root","read");
  TFile *fPt3 = new TFile("v2_D0method_pt6.5-7.5_y0.0-2.4_muPt0.0_centrality20-120.root","read");
  TFile *fPt4 = new TFile("v2_D0method_pt7.5-9.0_y0.0-2.4_muPt0.0_centrality20-120.root","read");
  TFile *fPt5 = new TFile("v2_D0method_pt9.0-12.0_y0.0-2.4_muPt0.0_centrality20-120.root","read");
  TFile *fPt6 = new TFile("v2_D0method_pt12.0-15.0_y0.0-2.4_muPt0.0_centrality20-120.root","read");
  TFile *fPt7 = new TFile("v2_D0method_pt15.0-50.0_y0.0-2.4_muPt0.0_centrality20-120.root","read");


  TH1D * NP1 = (TH1D*) fPt1 -> Get("hOutNP");
  TH1D * NP2 = (TH1D*) fPt2 -> Get("hOutNP");
  TH1D * NP3 = (TH1D*) fPt3 -> Get("hOutNP");
  TH1D * NP4 = (TH1D*) fPt4 -> Get("hOutNP");
  TH1D * NP5 = (TH1D*) fPt5 -> Get("hOutNP");
  TH1D * NP6 = (TH1D*) fPt6 -> Get("hOutNP");
  TH1D * NP7 = (TH1D*) fPt7 -> Get("hOutNP");

  double xNPpt1, yNPpt1, xNPpt2, yNPpt2, xNPpt3, yNPpt3, xNPpt4, yNPpt4, xNPpt5, yNPpt5, xNPpt6, yNPpt6, xNPpt7, yNPpt7;
  double xNPpt1Err, yNPpt1Err, xNPpt2Err, yNPpt2Err, xNPpt3Err, yNPpt3Err, xNPpt4Err, yNPpt4Err, xNPpt5Err, yNPpt5Err, xNPpt6Err, yNPpt6Err, xNPpt7Err, yNPpt7Err;

  yNPpt1= NP1->GetBinContent(0);
  yNPpt2= NP2->GetBinContent(0);
  yNPpt3= NP3->GetBinContent(0);
  yNPpt4= NP4->GetBinContent(0);
  yNPpt5= NP5->GetBinContent(0);
  yNPpt6= NP6->GetBinContent(0);
  yNPpt7= NP7->GetBinContent(0);
  yNPpt1Err= NP1->GetBinError(0);
  yNPpt2Err= NP2->GetBinError(0);
  yNPpt3Err= NP3->GetBinError(0);
  yNPpt4Err= NP4->GetBinError(0);
  yNPpt5Err= NP5->GetBinError(0);
  yNPpt6Err= NP6->GetBinError(0);
  yNPpt7Err= NP7->GetBinError(0);


  double xlow[] = {3.75, 5.5};
  double xlerr[] = {0.0, 0.0};
  
  double yprplow[] = {yNPpt1,yNPpt2}; // D0 method low pt
  double yprplowerr[] = {yNPpt1Err,yNPpt2Err}; // D0 method low pt

  double dyprplow[] = {0.077, 0.016}; // nominal prompt low pt
  double dyprplowerr[] = {0.028, 0.022}; // nominal prompt low pt

  double xbin[] = {7.0, 8.25, 10.5, 13.5, 32.5};
  double xerr[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double xerr2[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  double yprp[] = {yNPpt3,yNPpt4,yNPpt5,yNPpt6,yNPpt7}; // D0 method
  double yprperr[] = {yNPpt3Err,yNPpt4Err,yNPpt5Err,yNPpt6Err,yNPpt7Err}; // D0 method

  double dyprp[] = {0.030, 0.053, 0.051, 0.034, 0.021}; // nominal prompt 
  double dyprperr[] = {0.022, 0.015, 0.009, 0.010, 0.011}; // nominal prompt 

  TGraphErrors *v2ptPrp = new TGraphErrors(5, xbin, yprp, xerr, yprperr);
  TGraphErrors *v2ptdPrp = new TGraphErrors(5, xbin, dyprp, xerr, dyprperr);


  TGraphErrors *v2ptPrplow = new TGraphErrors(2, xlow, yprplow, xlerr, yprplowerr);
  TGraphErrors *v2ptdPrplow = new TGraphErrors(2, xlow, dyprplow, xlerr, dyprplowerr);

  TH1F *hPad = new TH1F("hPad",";p_{T} (GeV/c);v_{2}", 10, 0.0, 50.0);
  //hPad->SetMaximum(0.8);
  //hPad->SetMaximum(0.25);
  hPad->SetMaximum(0.2);
  //hPad->SetMinimum(-0.1);
  hPad->SetMinimum(-0.05);
  hPad->GetXaxis()->CenterTitle();
  hPad->GetYaxis()->CenterTitle();

  TCanvas *c1 = new TCanvas("c1","",880,800);
  c1->cd();
  hPad->Draw();

  v2ptdPrp->SetMarkerStyle(33);
  v2ptdPrp->SetMarkerColor(kViolet+7);
  v2ptdPrp->SetLineColor(kViolet+7);
  v2ptdPrp->SetMarkerSize(2.6);
  v2ptdPrp->Draw("p same");


  v2ptPrp->SetMarkerStyle(28);
  v2ptPrp->SetMarkerColor(kRed+1);
  v2ptPrp->SetLineColor(kRed+1);
  v2ptPrp->SetMarkerSize(2.6);
  v2ptPrp->Draw("p same");


  v2ptdPrplow->SetMarkerStyle(29);
  v2ptdPrplow->SetMarkerColor(kViolet+1);
  v2ptdPrplow->SetLineColor(kViolet+1);
  v2ptdPrplow->SetMarkerSize(2.6);
  v2ptdPrplow->Draw("p same");


  v2ptPrplow->SetMarkerStyle(46);
  v2ptPrplow->SetMarkerColor(kRed-4);
  v2ptPrplow->SetLineColor(kRed-4);
  v2ptPrplow->SetMarkerSize(2.6);
  v2ptPrplow->Draw("p same");

  TLatex *lt1 = new TLatex(); lt1->SetNDC();
  TLatex *CMS = new TLatex(); CMS->SetNDC();
  TLatex *Pre = new TLatex(); Pre->SetNDC();
  //CMS->SetTextSize(0.075*1.55);
  CMS->SetTextSize(0.04);
  CMS->SetTextFont(61);
  Pre->SetTextSize(0.04);
  Pre->SetTextFont(52);
  CMS->DrawLatex(0.60,0.87,"CMS");
  Pre->DrawLatex(0.70,0.87,"Preliminary");
	lt1->SetTextSize(0.040);
  lt1->DrawLatex(0.2,0.87,"Scalar Product");
  lt1->DrawLatex(0.62,0.94,"PbPb 1.7 nb^{-1} (5.02 TeV)");
  lt1->DrawLatex(0.2,0.74,"Cent. 10-60 \%");

  TLegend *leg1 = new TLegend(0.60,0.60,0.70,0.74);
  //TLegend *leg1 = new TLegend(0.55,0.66,0.68,0.85);
  //TLegend *leg1 = new TLegend(0.52,0.61,0.65,0.80);
  //TLegend *leg1 = new TLegend(0.52,0.55,0.65,0.74);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.03);

  //v2ptPrplowsys-> SetFillColorAlpha(kMagenta-10,0.3);
  //v2ptPrplowsys->SetFillStyle(3018);
  //v2ptPrpsys->SetFillStyle(3018);
  //v2ptPrpsys-> SetFillColorAlpha(kBlue-10,0.3);

  // Only prompt
  c1->cd();
  hPad->Draw();
  v2ptPrp->Draw("p same");
  v2ptdPrp->Draw("p same");
  v2ptPrplow->Draw("p same");
  v2ptdPrplow->Draw("p same");

  lt1->SetTextFont(42);
  lt1->DrawLatex(0.2,0.87,"Scalar Product");
  lt1->DrawLatex(0.2,0.80,"Cent. 10-60 \%");
  lt1->DrawLatex(0.2,0.73,"Non Prompt J/#psi");
  //lt1->DrawLatex(0.2,0.73,"|y| < 2.4");

  CMS_lumi_v2mass(c1,iPeriod,iPos);

  TLegend *leg2 = new TLegend(0.53,0.56,0.74,0.77);
  //TLegend *leg2 = new TLegend(0.60,0.57,0.79,0.77);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.03);
  //leg2->SetHeader("Prompt J/#psi");
  //leg2->AddEntry(v2ptPrp,"Prompt #psi(2S) |y| < 2.4");
  //leg2->AddEntry(v2ptdPrp,"Prompt #psi(2S) |y| < 2.4");
  leg2->AddEntry(v2ptdPrp,"Nominal method (|y|<2.4)");
  leg2->AddEntry(v2ptPrp,"D_{0} method (|y|<2.4)");
  leg2->AddEntry(v2ptdPrplow,"Nominal method (1.6<|y|<2.4)");
  leg2->AddEntry(v2ptPrplow,"D_{0} method (1.6<|y|<2.4)");
  leg2->Draw("same");

  c1->SaveAs("plot_D0_val_v2_nonprompt_jpsi_pt_weighted_7bin.png");
  c1->SaveAs("plot_D0_val_v2_nonprompt_jpsi_pt_weighted_7bin.pdf");
}
