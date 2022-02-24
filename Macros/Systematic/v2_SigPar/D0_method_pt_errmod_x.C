#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../../../cutsAndBin.h"
#include "../../../commonUtility.h"
#include "../../../JpsiUtility.h"

valErr getYield(float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0, int cLow=0, int cHigh=0, double v2=0);
double getFrac(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh, double v2);
double getFracErr(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh, double v2);
void D0_method_pt_errmod_x(
    float ptLow=3, float ptHigh=4.5,
    float yLow=1.6, float yHigh=2.4,
    int cLow=20, int cHigh=120
    )
{
  using namespace std;
  using namespace RooFit; 

  gStyle->SetOptStat(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.05);


  TString kineLabel = getKineLabel (ptLow, ptHigh,yLow, yHigh, 0.0, cLow, cHigh);
  TString DATE=Form("%i_%i",cLow/2,cHigh/2);

  //const int nV2Bins=18;
  //double v2Bin[nV2Bins+1] = {-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0.,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7};
  const int nV2Bins=16;
  double v2Bin[nV2Bins+1] = {-2.7,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0.,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.7};
  //const int nV2Bins=12;
  //double v2Bin[nV2Bins+1] = {-2.7,-1.8,-1.2,-0.9,-0.6,-0.3,0.,0.3,0.6,0.9,1.2,1.8,2.7};
  double frac[nV2Bins];
  double fracErr[nV2Bins];

  TFile* fNom = new TFile("/Users/hwan/tools/2019/CMS/JPsi/Jpsi_v2_PbPb2018/Macros/21_11_11/v2_D0method_pt3.0-4.5_y1.6-2.4_muPt0.0_centrality20-120_new.root");
  TH1D* hPRNom = (TH1D*)fNom->Get("hPR");

  TH1D* hPR = new TH1D("hPR",";#frac{Q_{2}Q_{2A}^{*}}{#sqrt{#frac{#LTQ_{2A}Q_{2B}^{*}#GT#LTQ_{2B}Q_{2C}^{*}#GT}{#LTQ_{2A}Q_{2C}^{*}#GT}}};",nV2Bins,v2Bin);
  TH1D* hNP = new TH1D("hNP",";#frac{Q_{2}Q_{2A}^{*}}{#sqrt{#frac{#LTQ_{2A}Q_{2B}^{*}#GT#LTQ_{2B}Q_{2C}^{*}#GT}{#LTQ_{2A}Q_{2C}^{*}#GT}}};",nV2Bins,v2Bin);
  TH1D* hPRbfr = new TH1D("hPRbfr",";#frac{Q_{2}Q_{2A}^{*}}{#sqrt{#frac{#LTQ_{2A}Q_{2B}^{*}#GT#LTQ_{2B}Q_{2C}^{*}#GT}{#LTQ_{2A}Q_{2C}^{*}#GT}}};b fraction",nV2Bins,v2Bin);
  TH1D* hNPbfr = new TH1D("hNPbfr",";#frac{Q_{2}Q_{2A}^{*}}{#sqrt{#frac{#LTQ_{2A}Q_{2B}^{*}#GT#LTQ_{2B}Q_{2C}^{*}#GT}{#LTQ_{2A}Q_{2C}^{*}#GT}}};b fraction",nV2Bins,v2Bin);
  TF1* f1;
  f1 = new TF1("f1","gaus",-2.7,2.7);
  for(int iv2=0;iv2<nV2Bins;iv2++) {
    frac[iv2]=getFrac(ptLow,ptHigh,yLow,yHigh,cLow,cHigh,v2Bin[iv2]);
    fracErr[iv2]=getFracErr(ptLow,ptHigh,yLow,yHigh,cLow,cHigh,v2Bin[iv2]);
    cout << v2Bin[iv2] << " - " << v2Bin[iv2+1] << " frac : " << frac[iv2] <<" +/- "<<fracErr[iv2]<< endl;
  }
  for(int iv2=0;iv2<nV2Bins;iv2++) {
    valErr yieldAA;
    yieldAA = getYield(ptLow,ptHigh,yLow,yHigh,cLow,cHigh,v2Bin[iv2]);
    cout << v2Bin[iv2] << " - " << v2Bin[iv2+1] << " yield : " << yieldAA.val <<" +/- "<< yieldAA.err << endl;
    double err1 = yieldAA.err;
    double err2 = fracErr[iv2];
    double errprop_npr = yieldAA.val*frac[iv2]*sqrt((err1/yieldAA.val)*(err1/yieldAA.val) + (err2/frac[iv2])*(err2/frac[iv2]));
    double errprop_prp = yieldAA.val*(1-frac[iv2])*sqrt((err1/yieldAA.val)*(err1/yieldAA.val) + (err2/(1-frac[iv2]))*(err2/(1-frac[iv2])));
    hPR->SetBinContent(iv2+1,yieldAA.val*(1-frac[iv2]));
    hPR->SetBinError(iv2+1,errprop_prp);
    hNP->SetBinContent(iv2+1,yieldAA.val*(frac[iv2]));
    hNP->SetBinError(iv2+1,errprop_npr);
    hPRbfr->SetBinContent(iv2+1,frac[iv2]);
    hNPbfr->SetBinContent(iv2+1,frac[iv2]);
    hPRbfr->SetBinError(iv2+1,fracErr[iv2]);
    hNPbfr->SetBinError(iv2+1,fracErr[iv2]);
    if(iv2==0 || iv2==nV2Bins-1) {
      hPR->SetBinContent(iv2+1,yieldAA.val*(1-frac[iv2])/2.0);
      hPR->SetBinError(iv2+1,errprop_prp/2.0);
      hNP->SetBinContent(iv2+1,yieldAA.val*(frac[iv2])/2.0);
      hNP->SetBinError(iv2+1,errprop_npr/2.0);
      //        if(iv2==1 || iv2==nV2Bins-2) {
      //            hPR->SetBinContent(iv2+1,yieldAA.val*(1-frac[iv2])/2.0);
      //            hPR->SetBinError(iv2+1,errprop_prp/2.0);
      //            hNP->SetBinContent(iv2+1,yieldAA.val*(frac[iv2])/2.0);
      //            hNP->SetBinError(iv2+1,errprop_npr/2.0);
    }
    }
  cout << "Total : " << hPR->Integral() << endl;
  cout << "Mean : " << hPR->GetMean() << endl;
  cout << "Mean Error : " << hPR->GetMeanError() << endl;
  int text_size=17;
  double text_x = 0.18;
  double text_y = 0.8;
  TCanvas *c1 = new TCanvas("c1",Form("%.1f <p_{T} < %.1f (GeV/c)",ptLow,ptHigh),1440,700);
  c1->Divide(2,1);
  c1->cd(1);
  c1->cd(1)->SetBottomMargin(0.2);
  hPR->SetTitle(Form("%.1f < p_{T} < %.1f (GeV/c)",ptLow,ptHigh));
  hPR->GetXaxis()->CenterTitle();
  hPR->GetXaxis()->SetTitleOffset(2.2);
  hPR->GetXaxis()->SetTitleSize(0.03);
  hPR->Draw("PE");
  drawText("Prompt J/#psi",text_x,text_y,text_color,text_size+2);
  drawText(Form("Entries = %.f",hPR->GetEntries()),text_x,text_y-0.05,text_color,text_size);
  drawText(Form("Mean = %.4f #pm %.4f",hPR->GetMean(),hPR->GetMeanError()),text_x,text_y-0.1,text_color,text_size);
  c1->cd(2);
  c1->cd(2)->SetBottomMargin(0.2);
  hNP->GetXaxis()->CenterTitle();
  hNP->GetXaxis()->SetTitleOffset(2.2);
  hNP->GetXaxis()->SetTitleSize(0.03);
  hNP->Draw("PE");
  drawText("NonPrompt J/#psi",text_x,text_y,text_color,text_size+2);
  drawText(Form("Entries = %d",nV2Bins),text_x,text_y-0.05,text_color,text_size);
  drawText(Form("Mean = %.4f #pm %.4f",hNP->GetMean(),hNP->GetMeanError()),text_x,text_y-0.1,text_color,text_size);
  c1->SaveAs(Form("./SigPar_x/plot_D0method_x_%s.png",kineLabel.Data()));
  c1->SaveAs(Form("./SigPar_x/plot_D0method_x_%s.pdf",kineLabel.Data()));


  TF1 *line1 = new TF1("line1","[0]",-2.8,2.8);
  line1->SetLineColor(kRed+2);
  line1->SetLineStyle(2);

  TCanvas *c2 = new TCanvas("c2",",",660,600);
  c2->cd();
  c2->SetBottomMargin(0.2);
  hPRbfr->Fit(line1,"rm");
  //if(ptLow==3.0) hPRbfr->GetYaxis()->SetRangeUser(0.0,0.3);
  //else if(ptLow==9.0) hPRbfr->GetYaxis()->SetRangeUser(0.2,0.5);
  //else if(ptLow==12.0) hPRbfr->GetYaxis()->SetRangeUser(0.3,0.6);
  //else if(ptLow==15.0) hPRbfr->GetYaxis()->SetRangeUser(0.3,0.7);
  //else hPRbfr->GetYaxis()->SetRangeUser(0.1,0.4);
  double mean = line1->GetParameter(0);
  hPRbfr->GetYaxis()->SetRangeUser(0.0,1.0);
  //hPRbfr->GetXaxis()->SetTitle("v_{2}");
  hPRbfr->GetXaxis()->CenterTitle();
  hPRbfr->GetXaxis()->SetTitleOffset(2.2);
  hPRbfr->GetXaxis()->SetTitleSize(0.03);
  hPRbfr->Draw("PE");
  //drawText("Prompt J/#psi",text_x,text_y,text_color,text_size+2);
  //drawText(Form("Entries = %.f",hPR->GetEntries()),text_x,text_y-0.05,text_color,text_size);
  //drawText(Form("Mean = %.4f #pm %.4f",hPR->GetMean(),hPR->GetMeanError()),text_x,text_y-0.1,text_color,text_size);
  drawText(Form("%0.1f < p_{T} < %0.1f (GeV/c)",ptLow,ptHigh),text_x,text_y+0.3,text_color,text_size);
  drawText(Form("mean b fraction : %0.3f #pm %0.3f ",line1->GetParameter(0),line1->GetParError(0)),text_x,text_y-0.1,text_color,text_size);

  c2->SaveAs(Form("./SigPar_x/plot_D0method_bfractions_x_%s.png",kineLabel.Data()));
  c2->SaveAs(Form("./SigPar_x/plot_D0method_bfractions_x_%s.pdf",kineLabel.Data()));

  TCanvas *c3 = new TCanvas("c3","",660,600);
  c3->cd();
  hPRNom->Divide(hPR);
  for(int ibin =1; ibin<= hPRNom->GetNbinsX(); ibin++){
	  hPRNom->SetBinError(ibin, 0); }
  hPRNom->GetYaxis()->SetTitle("#frac{Nom}{x}");
  hPRNom->GetYaxis()->CenterTitle();
  hPRNom->GetYaxis()->SetRangeUser(0.7,1.3);
  hPRNom->Draw("P");
  jumSun(-2.8,1,2.8,1);
  c3->SaveAs("./Ratio_x.png");

  TH1D *hOutPR = new TH1D("hOutPR","",1,0,1);
  hOutPR->SetBinContent(0,hPR->GetMean());
  hOutPR->SetBinError(0,hPR->GetMeanError());

  TH1D *hOutNP = new TH1D("hOutNP","",1,0,1);
  hOutNP->SetBinContent(0,hNP->GetMean());
  hOutNP->SetBinError(0,hNP->GetMeanError());

  TFile *outFile;
  outFile = new TFile(Form("./SigPar_x/v2_D0method_%s_x.root",kineLabel.Data()),"recreate");
  outFile->cd();
  hOutPR->Write();
  hOutNP->Write();

}

//Get Yield
valErr getYield(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh, double v2) {
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);
  TFile* inf = new TFile(Form("./roots/2DFit_10_60/Mass_v2Bins_x/Mass_FixedFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1_v2_%.2f.root", kineLabel.Data(),v2));
  //TFile* inf = new TFile(Form("./roots/2DFit_0_180/Mass_v2Bins_x/Mass_FixedFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1_v2_%.2f.root", kineLabel.Data(),v2));
  //TFile* inf = new TFile(Form("../Macros/2021_04_22/roots/2DFit_210915/Mass/MassFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel.Data()));
  //RooWorkspace* ws = (RooWorkspace*)inf->Get("workspace");
  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  valErr ret;
  ret.val = fitResults->GetBinContent(1);
  ret.err = fitResults->GetBinError(1);
  //cout << Form("./roots/2DFit_10_60/Mass/Mass_constraintFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1_v2_%.2f.root", kineLabel.Data(),v2) << endl;
  //cout << Form("../Macros/2021_04_22/roots/2DFit_210915/Mass/MassFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel.Data()) << endl;
  //cout << kineLabel << ": " << " & " << ret.val << " $\pm$ " << ret.err << " & " <<ws->var("nBkg")->getVal() << " $\pm$ "<< ws->var("nBkg")->getError() << "\\\\" << endl;
  return ret;
}
double getFrac(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh, double v2) {
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);
  //TFile* inf = new TFile(Form("../roots/2DFit_210915/Final/2DFitResult_Inclusive_%s.root", kineLabel.Data()));
  //TFile* inf = new TFile(Form("./roots/2DFit_10_60/Final_v2Bins_x/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1_v2_%.2f.root", kineLabel.Data(),v2));
  TFile* inf = new TFile(Form("../../21_11_11/roots/2DFit_10_60/Final_v2Bins/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1_v2_%.2f.root", kineLabel.Data(),v2));
  //TFile* inf = new TFile(Form("./roots/2DFit_0_180/Final_v2Bins_x/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1_v2_%.2f.root", kineLabel.Data(),v2));
  TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");
  double frac;
  frac = fitResults->GetBinContent(1);
  return frac;
}
double getFracErr(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh, double v2) {
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);
  //TFile* inf = new TFile(Form("../roots/2DFit_210915/Final/2DFitResult_Inclusive_%s.root", kineLabel.Data()));
  //TFile* inf = new TFile(Form("./roots/2DFit_10_60/Final_v2Bins_x/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1_v2_%.2f.root", kineLabel.Data(),v2));
  TFile* inf = new TFile(Form("../../21_11_11/roots/2DFit_10_60/Final_v2Bins/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1_v2_%.2f.root", kineLabel.Data(),v2));
  //TFile* inf = new TFile(Form("./roots/2DFit_0_180/Final_v2Bins_x/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1_v2_%.2f.root", kineLabel.Data(),v2));
  TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");
  double frac;
  frac = fitResults->GetBinError(1);
  return frac;
}
