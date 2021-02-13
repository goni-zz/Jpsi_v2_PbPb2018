void drawsame() {
	TFile * f1 = new TFile("mc_eff_vs_pt_cent_prompt_pbpb_Jpsi_fit_tnp.root","read");
	//TH1D * h1=f1.Get("hpt_eff_Trig_1");
	TObject * h1=f1->Get("hpt_eff_Trig_1");
	TObject * h2=f1->Get("hpt_eff_Trig_2");
	TObject * h3=f1->Get("hpt_eff_Trig_3");
	TFile * f2 = new TFile("mc_eff_vs_pt_cent_nprompt_pbpb_Jpsi_fit_tnp_2.root","read");
	//TH1D * h2=f2.Get("hpt_eff_Trig_1");
	TObject * h6=f2->Get("hpt_eff_Trig_6");
	TObject * h7=f2->Get("hpt_eff_Trig_7");
	TObject * h8=f2->Get("hpt_eff_Trig_8");

	TLatex *lt1 = new TLatex();
	lt1->SetNDC();
	lt1->SetTextSize(0.03);
	auto legend = new TLegend(0.6,0.84);
	auto legend2 = new TLegend(0.6,0.84);

	gStyle->SetOptFit(0);
	TCanvas * cpt_eff1 = new TCanvas("cpt_eff1","cpt_eff1",0,0,900,800);
	cpt_eff1->cd();
	h1->Draw("E");
	h6->Draw("same");
	legend->AddEntry("hpt_eff_Trig_1","|y|: 0.0-2.4, 0-90%, 6.5<Pt<50 Prompt","lep");
	legend->AddEntry("hpt_eff_Trig_6","|y|: 0.0-2.4, 0-90%, 6.5<Pt<50 NPrompt","lep");
	legend->SetBorderSize(0);
	legend->Draw("same");
	lt1->SetTextSize(0.03);

	TCanvas * cpt_eff2 = new TCanvas("cpt_eff2","cpt_eff2",0,0,900,800);
	cpt_eff2->cd();
	h2->Draw("E");
	h3->Draw("same");
	h7->Draw("same");
	h8->Draw("same");
	legend2->AddEntry("hpt_eff_Trig_2","|y|: 0.0-1.6, 0-90%, 6.5<Pt<50 Prompt","lep");
	legend2->AddEntry("hpt_eff_Trig_3","|y|: 1.6-2.4, 0-90%, 3<Pt<50 Prompt","lep");
	legend2->AddEntry("hpt_eff_Trig_7","|y|: 0.0-1.6, 0-90%, 6.5<Pt<50 NPrompt","lep");
	legend2->AddEntry("hpt_eff_Trig_8","|y|: 1.6-2.4, 0-90%, 3<Pt<50 NPrompt","lep");
	legend2->SetBorderSize(0);
	legend2->Draw("same");
	lt1->SetTextSize(0.03);

	cpt_eff1->Write();
	cpt_eff2->Write();

	TFile* outFile = new TFile("Drawsame.root","RECREATE");
}
