#include <iostream>
#include "../Style.h"
#include "../commonUtility.h"
#include "../cutsAndBin.h"
#include "../CMS_lumi_square.C"
#include "../tdrstyle.C"

void draw_V2_pt(){
	setTDRStyle();
	writeExtraText = true;
	int iPeriod = 2;
	int iPos = 33;

	const int nBin = 5;
	
	double ptBin[nBin+1] = {3,4.5,6.5,9,12,30};

	TH1D *hist_v2;
	hist_v2 = new TH1D("hist_v2",";p_{T}^{J#Psi} (GeV/C);v_{2}",nBin,ptBin);

	hist_v2->SetBinContent(1,0.056);
	hist_v2->SetBinContent(2,0.038);
	hist_v2->SetBinContent(3,0.056);
	hist_v2->SetBinContent(4,0.047);
	hist_v2->SetBinContent(5,0.034);
	
	hist_v2->SetBinError(1,0.020);
	hist_v2->SetBinError(2,0.012);
	hist_v2->SetBinError(3,0.011);
	hist_v2->SetBinError(4,0.014);
	hist_v2->SetBinError(5,0.016);

	hist_v2->GetXaxis()->CenterTitle();
	hist_v2->GetYaxis()->CenterTitle();

	TCanvas *c1;
    TLegend *leg;
    TLatex* globtex = new TLatex();
    globtex->SetNDC();
    globtex->SetTextAlign(12); // left-center
    globtex->SetTextFont(42);
    globtex->SetTextSize(0.040);
    double sz_init = 0.867; double sz_step = 0.0535;
    c1=new TCanvas("c1","c1",600,600);
    c1->cd();
    hist_v2->SetMarkerStyle(kFullCircle);
	hist_v2->SetMarkerColor(kBlue+2);
	hist_v2->SetLineColor(kBlue+2);
    hist_v2->GetXaxis()->CenterTitle();
    hist_v2->GetYaxis()->CenterTitle();
    hist_v2->GetYaxis()->SetTitleOffset(1.2);
    hist_v2->GetXaxis()->SetTitleOffset(1.);
    hist_v2->GetYaxis()->SetTitle("v_{2}^{J/#psi}");
    hist_v2->GetXaxis()->SetTitle("p_{T}^{J/#psi} (GeV/c)");
    hist_v2->SetMinimum(-0.05);
    hist_v2->SetMaximum(0.15);
    hist_v2->GetXaxis()->SetRangeUser(0,30);
    hist_v2->GetYaxis()->SetRangeUser(-0.05,0.15);
    hist_v2->Draw("PE");

    TF1 *line = new TF1("line","0",-0.5,30);
    line->SetLineWidth(1);
    line->SetLineStyle(7);
    line->SetLineColor(kBlack);

    line->Draw("same");


    TString perc = "%";
    leg = new TLegend(0.74, 0.70, 0.925, 0.77);
    SetLegendStyle(leg);
    leg->SetTextSize(0.044);
    leg->AddEntry(hist_v2,"J/#psi","pe");
    leg->Draw("same");
    globtex->DrawLatex(0.23, sz_init, "1.6 < |y| < 2.4");
    globtex->DrawLatex(0.23, sz_init-sz_step, Form("Cent. 10-60%s",perc.Data()));
    CMS_lumi_square(c1,iPeriod,iPos);
    c1->Update();
    c1->SaveAs("v2_pt_Cent1060.pdf");


}
