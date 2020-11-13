#include <iostream>
#include "../rootFitHeaders.h"
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBin.h"
#include "../CMS_lumi_v2mass.C"
#include "../tdrstyle.C"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit;

void MassFit(
		float ptLow=3, float ptHigh=6.5,
		float yLow=1.6, float yHigh=2.4,
		int cLow=0, int cHigh=200,
		float muPtCut=0.0,
		bool whichModel=0,  // Nominal=0, Alternative=1
		int ICset=1
		)
{
	gStyle->SetEndErrorSize(0);
	gSystem->mkdir("../roots/2DFit/");
	gSystem->mkdir("../figs/2DFit/");

    RooMsgService::instance().getStream(0).removeTopic(Caching);
    RooMsgService::instance().getStream(1).removeTopic(Caching);
    RooMsgService::instance().getStream(0).removeTopic(Plotting);
    RooMsgService::instance().getStream(1).removeTopic(Plotting);
    RooMsgService::instance().getStream(0).removeTopic(Integration);
    RooMsgService::instance().getStream(1).removeTopic(Integration);
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

	TFile* f1; TFile* f2; TFile* f3;
	TString kineCut;
	TString SigCut;
	TString BkgCut;

	f1 = new TFile("../skimmedFiles/OniaRooDataSet_isMC0_JPsi1SW_20200928.root");
	kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f",ptLow, ptHigh, yLow, yHigh);
	SigCut  = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.8 && mass<3.2",ptLow, ptHigh, yLow, yHigh);
	BkgCut  = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && ((mass>2.6 && mass <= 2.8) || (mass>=3.2&&mass<3.5))",ptLow, ptHigh, yLow, yHigh);
	TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

	kineCut = accCut+kineCut;
	SigCut = SigCut;
	BkgCut = BkgCut;

	RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");

	RooWorkspace *ws = new RooWorkspace("workspace");
	ws->import(*dataset);
	cout << "####################################" << endl;
	RooDataSet *reducedDS_A = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), SigCut.Data() );
	RooDataSet *reducedDS_B = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), BkgCut.Data() );
	reducedDS_A->SetName("reducedDS_A");
	reducedDS_B->SetName("reducedDS_B");

	RooCategory tp("tp","tp");
	tp.defineType("A");
	tp.defineType("B");

	RooDataSet* dsAB = new RooDataSet("dsAB","dsAB",RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))),Index(tp),Import("A",*reducedDS_A),Import("B",*reducedDS_B));
	cout << "******** New Combined Dataset ***********" << endl;
	dsAB->Print();
	ws->import(*dsAB);
	ws->var("mass")->setRange(massLow, massHigh);
	ws->var("mass")->Print();


//***********************************************************************
//****************************** MASS FIT *******************************
//***********************************************************************

	//         The order is {sigma_1,  x, alpha_1, n_1,   f, err_mu, err_sigma, m_lambda}
	double paramsupper[8] = {0.2,    3.0,   5.321, 5.0, 1.0,   25.0,      25.0,     25.0};
	double paramslower[8] = {0.02,   0.0,     1.0, 1.0, 0.0,    0.0,       0.0,      0.0};
	//SIGNAL: initial params
	double sigma_1_init = 0.1;
	double x_init = 0.9;
	double alpha_1_init = 1.5;
	double n_1_init = 2.6;
	double f_init = 0.5;
	if (ICset>1 && ICset<4) {
		sigma_1_init = 0.3;
		x_init = 0.3;
		alpha_1_init = 2.6;
		n_1_init = 3.0;
		f_init = 0.1;
	}

	//SIGNAL
	RooRealVar    mean("m_{J/#Psi}","mean of the signal gaussian mass PDF",pdgMass.JPsi, pdgMass.JPsi -0.1, pdgMass.JPsi + 0.1 ) ;
	RooRealVar    sigma_1_A("sigma_1_A","width/sigma of the signal gaussian mass PDF",sigma_1_init, paramslower[0], paramsupper[0]);
	RooRealVar   *x_A = new RooRealVar("x_A","sigma ratio ", x_init, paramslower[1], paramsupper[1]);
	RooFormulaVar sigma_2_A("sigma_2_A","@0*@1",RooArgList(sigma_1_A, *x_A) );
	RooRealVar    alpha_1_A("alpha_1_A","tail shift", alpha_1_init, paramslower[2], paramsupper[2]);
	RooFormulaVar alpha_2_A("alpha_2_A","1.0*@0",RooArgList(alpha_1_A) );
	RooRealVar    n_1_A("n_1_A","power order", n_1_init , paramslower[3], paramsupper[3]);
	RooFormulaVar n_2_A("n_2_A","1.0*@0",RooArgList(n_1_A) );
	RooRealVar   *f = new RooRealVar("f","cb fraction", f_init, paramslower[4], paramsupper[4]);
	//Set up crystal ball shapes
	RooCBShape* cb_1_A = new RooCBShape("cball_1_A", "cystal Ball", *(ws->var("mass")), mean, sigma_1_A, alpha_1_A, n_1_A);
	RooAddPdf*  pdfMASS_Jpsi;
	//DOUBLE CRYSTAL BALL
	RooCBShape* cb_2_A = new RooCBShape("cball_2_A", "cystal Ball", *(ws->var("mass")), mean, sigma_2_A, alpha_2_A, n_2_A);
	pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );
	RooRealVar *N_Jpsi= new RooRealVar("N_Jpsi","inclusive Jpsi signals",0,1000000);
	pdfMASS_Jpsi = new RooAddPdf("pdfMASS_Jpsi","Signal ",RooArgList(*cb_1_A,*cb_2_A), RooArgList(*f) );
	//BACKGROUND
	double err_sigma_init = 5;
	double err_mu_init = 8;
	double m_lambda_init = 5;
	if (ICset>2) {
		err_mu_init = 5;
		m_lambda_init = 5;
	}

	//A: Mass Bkg function
	RooRealVar err_mu_A("#mu_A","err_mu", err_mu_init,  paramslower[5], paramsupper[5]) ;
	RooRealVar err_sigma_A("#sigma_A","err_sigma", err_sigma_init, paramslower[6], paramsupper[6]);
	RooRealVar m_lambda_A("#lambda_A","m_lambda",  m_lambda_init, paramslower[7], paramsupper[7]);
	//THIS IS THE BACKGROUND FUNCTION
	//RooGenericPdf *pdfMASS_bkg = new RooGenericPdf("pdfMASS_bkg","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_A));
	//RooRealVar *N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",10000,0,1000000);

	RooRealVar *sl1 = new RooRealVar("sl1","sl1",0.6,0.,1);
	RooRealVar *cnst1 = new RooRealVar("cnst1","cnst1",0.8,0.,1);
	//RooGenericPdf *bkg = new RooGenericPdf("bkg","Background","@0*@1+@2",RooArgList( *(ws->var("mass")), sl1, cnst1) );
	RooChebychev *pdfMASS_bkg;
	pdfMASS_bkg = new RooChebychev("pdfMASS_bkg","Background",*(ws->var("mass")),RooArgList(*sl1,*cnst1));
	RooRealVar *N_Bkg = new RooRealVar("N_Bkg","fraction of component 1 in bkg",10000,0,1000000);
	//Build the model
	//Model A: Mass
	RooAddPdf* pdfMASS_Tot = new RooAddPdf();
	pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *pdfMASS_bkg),RooArgList(*N_Jpsi,*N_Bkg));
	//pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","Jpsi + Bkg",RooArgList(*pdfMASS_Jpsi, *bkg_1order),RooArgList(*N_Jpsi,*N_Bkg));
	//pdfMASS_Tot = new RooAddPdf("pdfMASS_Tot","PR Jpsi + NP Jpsi + Bkg",RooArgList(*cb_1_A, *cb_2_A, *bkg),RooArgList(*N_JpsiPR,*N_JpsiNP,*N_Bkg));
	ws->import(*pdfMASS_Tot);
	cout << endl << "********* Starting Mass Dist. Fit **************" << endl << endl;
	RooFitResult* fitMass = ws->pdf("pdfMASS_Tot")->fitTo(*dsAB,Save(), Hesse(kTRUE), Range(massLow,massHigh), Timer(kTRUE), Extended(kTRUE), NumCPU(nCPU), PrintLevel(-1));
	cout << endl << "********* Finished Mass Dist. Fit **************" << endl << endl;

	TCanvas* c_A =  new TCanvas("canvas_A","My plots",4,4,550,520);
	c_A->cd();
	TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0, 0.16, 0.98, 1.0);
	pad_A_1->SetTicks(1,1);
	pad_A_1->Draw(); pad_A_1->cd();
	RooPlot* myPlot_A = ws->var("mass")->frame(nMassBin); // bins
	myPlot_A->SetTitle("");
	ws->data("dsAB")->plotOn(myPlot_A,Name("dataHist_A"));

	pad_A_1->cd();
	gPad->SetLogy();
	RooPlot* myPlot2_A = (RooPlot*)myPlot_A->Clone();
	dsAB->plotOn(myPlot2_A,Name("dataOS"),MarkerSize(.8));
	ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_tot"), LineColor(kBlack));
	//ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("Sig_A"),Components(RooArgSet(*pdfMASS_Jpsi)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
	ws->pdf("pdfMASS_Tot")->plotOn(myPlot2_A,Name("pdfMASS_Bkg"),Components(RooArgSet(*pdfMASS_bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
	//make a pretty plot
	myPlot2_A->SetFillStyle(4000);
	myPlot2_A->GetYaxis()->SetTitleOffset(1.43);
	//myPlot2_A->GetYaxis()->CenterTitle();
	//myPlot2_A->GetYaxis()->SetTitleSize(0.058);
	//myPlot2_A->GetYaxis()->SetLabelSize(0.054);
	myPlot2_A->GetYaxis()->SetRangeUser(ws->var("N_Jpsi")->getVal()/100, ws->var("N_Jpsi")->getVal());
	//myPlot2_A->SetMinimum(2*10);
	myPlot2_A->GetXaxis()->SetLabelSize(0);
	myPlot2_A->GetXaxis()->SetTitleSize(0);
	myPlot2_A->GetXaxis()->CenterTitle();
	myPlot2_A->GetXaxis()->SetRangeUser(massLow,massHigh);
	myPlot2_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
	myPlot2_A->Draw();

	drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
	if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
	else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
	drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);
	drawText(Form("n_{J/#psi} = %.f #pm %.f",ws->var("N_Jpsi")->getVal(),ws->var("N_Jpsi")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
	drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);

	TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0, 0.006, 0.98, 0.227);
	c_A->cd();
	pad_A_2->Draw();
	pad_A_2->cd();
	pad_A_2->SetTopMargin(0); // Upper and lower plot are joined
	pad_A_2->SetBottomMargin(0.67);
	pad_A_2->SetBottomMargin(0.4);
	pad_A_2->SetFillStyle(4000);
	pad_A_2->SetFrameFillStyle(4000);
	pad_A_2->SetTicks(1,1);

	RooPlot* frameTMP = (RooPlot*)myPlot2_A->Clone("TMP");
	RooHist* hpull_A = frameTMP->pullHist("dataOS","pdfMASS_tot", true);
	hpull_A->SetMarkerSize(0.8);
	RooPlot* pullFrame_A = ws->var("mass")->frame(Title("Pull Distribution")) ;
	pullFrame_A->addPlotable(hpull_A,"P") ;
	pullFrame_A->SetTitle("");
	pullFrame_A->SetTitleSize(0);
	pullFrame_A->GetYaxis()->SetTitleOffset(0.3) ;
	pullFrame_A->GetYaxis()->SetTitle("Pull") ;
	pullFrame_A->GetYaxis()->SetTitleSize(0.15) ;
	pullFrame_A->GetYaxis()->SetLabelSize(0.15) ;
	pullFrame_A->GetYaxis()->SetRangeUser(-3.8,3.8);
	pullFrame_A->GetYaxis()->CenterTitle();

	pullFrame_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
	pullFrame_A->GetXaxis()->SetTitleOffset(1.05) ;
	pullFrame_A->GetXaxis()->SetLabelOffset(0.04) ;
	pullFrame_A->GetXaxis()->SetLabelSize(0.15) ;
	pullFrame_A->GetXaxis()->SetTitleSize(0.15) ;
	pullFrame_A->GetXaxis()->CenterTitle();

	pullFrame_A->GetYaxis()->SetTickSize(0.04);
	pullFrame_A->GetYaxis()->SetNdivisions(404);
	pullFrame_A->GetXaxis()->SetTickSize(0.03);
	pullFrame_A->Draw() ;

	TLine *l1 = new TLine(massLow,0,massHigh,0);
	l1->SetLineStyle(1);
	l1->Draw("same");

	printChi2(ws, pad_A_2, frameTMP, fitMass, "mass", "dataOS", "pdfMASS_tot", nMassBin, false);

	RooRealVar *Par1 = ws->var("sl1");
	RooRealVar *Par2 = ws->var("cnst1");
	cout << "Chebychev Par1 : " << Par1->getVal() << " +/- " << Par1->getError() << endl;
	cout << "Chebychev Par2 : " << Par2->getVal() << " +/- " << Par2->getError() << endl;

	ws->Print();
	RooArgSet* fitargs = new RooArgSet();
	fitargs->add(fitMass->floatParsFinal());
	//RooDataSet *datasetMass = new RooDataSet("datasetMass","dataset with Mass Fit result", RooArgList(*ws->var("N_Jpsi"),*ws->var("N_Bkg")) );
	RooDataSet *datasetMass = new RooDataSet("datasetMass","dataset with Mass Fit result", *fitargs );
	datasetMass->add(*fitargs);
	RooWorkspace *wsmass = new RooWorkspace("workspaceMass");

	c_A->Update();
	TString kineLabel = getKineLabel (ptLow, ptHigh,yLow, yHigh, muPtCut, cLow, cHigh);
	c_A->SaveAs(Form("../figs/2DFit/Mass_%s.pdf",kineLabel.Data()));

	TFile* outFile;
	outFile = new TFile(Form("../roots/2DFit/MassFitResult_%s.root",kineLabel.Data()),"recreate");
	pdfMASS_Tot->Write();
	datasetMass->Write();
	outFile->Close();

	fitMass->Print();
}

