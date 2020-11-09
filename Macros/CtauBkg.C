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

void CtauBkg(
		float ptLow=3, float ptHigh=6.5,
		float yLow=1.6, float yHigh=2.4,
		int cLow=0, int cHigh=200,
		float muPtCut=0.0,
		bool whichModel=0,  // Nominal=0, Alternative=1
		int ICset=1
		)
{
	gStyle->SetEndErrorSize(0);
	TFile* f1; TFile* fMass; TFile* fCErr; TFile* fCRes;
	TString kineCut;
	TString SigCut;
	TString BkgCut;
	TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

	f1 = new TFile("../skimmedFiles/OniaRooDataSet_isMC0_JPsi1SW_2020819.root");
	fMass = new TFile(Form("../roots/2DFit/MassFitResult_%s.root",kineLabel.Data()));
	fCErr = new TFile(Form("../roots/2DFit/CtauErrResult_%s.root",kineLabel.Data()));
	fCRes = new TFile(Form("../roots/2DFit/CtauResResult_%s.root",kineLabel.Data()));

	kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f",ptLow, ptHigh, yLow, yHigh);
	SigCut  = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.8 && mass<3.2",ptLow, ptHigh, yLow, yHigh);
	BkgCut  = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && ((mass>2.6 && mass <= 2.8) || (mass>=3.2&&mass<3.5))",ptLow, ptHigh, yLow, yHigh);
	TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

	kineCut = accCut+kineCut;
	SigCut = SigCut;
	BkgCut = BkgCut;


	RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
	RooDataSet *datasetMass = (RooDataSet*)fMass->Get("datasetMass");
	RooAddPdf* pdfMASS_Tot = (RooAddPdf*)fMass->Get("pdfMASS_Tot");
	RooDataSet *dataw_Bkg = (RooDataSet*)fCErr->Get("dataw_Bkg");
	RooDataSet *dataw_Sig = (RooDataSet*)fCErr->Get("dataw_Sig");
	RooAddPdf* GaussModel_Tot = (RooAddPdf*)fCRes->Get("GaussModel_Tot");

	RooWorkspace *ws = new RooWorkspace("workspace");
	ws->import(*dataset);
	ws->import(*datasetMass);
	ws->import(*pdfMASS_Tot);
	ws->import(*dataw_Bkg);
	ws->import(*dataw_Sig);
	ws->import(*GaussModel_Tot);
	cout << "####################################" << endl;
	RooDataSet *reducedDS_A = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("N_Jpsi"))), SigCut.Data() );
	RooDataSet *reducedDS_B = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("N_Bkg"))), BkgCut.Data() );
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
	ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
	ws->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);
	ws->var("ctau3DErr")->setRange(ctauErrLow, ctauErrHigh);
	ws->var("ctau3DErr")->setRange("ctauErrRange",ctauErrLow, ctauErrHigh);
	ws->var("ctau3DRes")->setRange(ctauResLow, ctauResHigh);
	ws->var("ctau3DRes")->setRange("ctauResRange", ctauResLow, ctauResHigh);
	ws->var("mass")->Print();
	ws->var("ctau3D")->Print();
	ws->var("ctau3DErr")->Print();

	//***********************************************************************
	//**************************** Bkg CTAU FIT *****************************
	//***********************************************************************
	cout << endl << "************** Start BKG Ctau Fit *****************" << endl << endl;
	//make parameter 3 exp
	ws->factory("zeroMean[0.]");
	ws->factory("fDFSS[0.5, 1e-6, 1.1]");
	ws->factory("fDLIV[0.5, 1e-6, 1.]");
	ws->factory("lambdaDDS_Bkg[0.9, 1e-6, 1.1]");
	ws->factory("lambdaDF_Bkg[0.3, 1e-6, 1.1]");
	ws->factory("lambdaDSS_Bkg[0.8, 1e-6, 1.1]");
	//parameters fixed by Resolution model
	ws->var("ctauRes_mean")->setConstant(kTRUE);
	ws->var("s1_CtauRes")->setConstant(kTRUE);
	ws->var("s2_CtauRes")->setConstant(kTRUE);
	ws->var("s3_CtauRes")->setConstant(kTRUE);
	ws->var("f2_CtauRes")->setConstant(kTRUE);
	ws->var("f_CtauRes")->setConstant(kTRUE);
	ws->var("zeroMean")->setConstant(kTRUE);
	//make res model
	ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes1", "ctau3D",
				"ctauRes_mean", //"ctau1_CtauRes",
				"s1_CtauRes",
				"zeroMean",
				"ctau3DErr"
				));
	ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes2", "ctau3D",
				"ctauRes_mean", //"ctau2_CtauRes",
				"s2_CtauRes",
				"zeroMean",
				"ctau3DErr"
				));
	ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes3", "ctau3D",
				"ctauRes_mean", //"ctau3_CtauRes",
				"s3_CtauRes",
				"zeroMean",
				"ctau3DErr"
				));
	ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "ctauRes32", "ctauRes3", "ctauRes2", "f2_CtauRes"));
	ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "pdfCTAURES", "ctauRes32", "ctauRes1", "f_CtauRes"));
	//make 3 exp
	ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS", "ctau3D",
				"lambdaDSS_Bkg", "pdfCTAURES"));
	ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF", "ctau3D",
				"lambdaDF_Bkg", "pdfCTAURES"));
	ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", "pdfCTAUDDS", "ctau3D",
				"lambdaDDS_Bkg", "pdfCTAURES"));
	ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU1",
				"fDFSS", "pdfCTAUDSS", "pdfCTAUDF"));
	ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_BkgNoPR",//NP
				"fDLIV", "pdfCTAU1", "pdfCTAUDDS"));
	ws->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_BkgPR",//PR
				"pdfCTAURES"//resolution model
				));
	//b_Bkg fraction
	ws->factory("b_Bkg[0.6, 0., 1.]");//NP fraction for bkg
	ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_Bkg",
				"b_Bkg",
				"pdfCTAUCOND_BkgNoPR",
				"pdfCTAUCOND_BkgPR"
				));
	ws->factory(Form("RooExtendPdf::%s(%s,%s)", "pdfTot_Bkg", "pdfCTAUCOND_Bkg", "N_Bkg"));//N_Bkg is number of bkg from dataw_Bkg
	RooAddPdf* pdfTot_Bkg = new RooAddPdf("pdfTot_Bkg");
	//RooAbsPdf* ctauBkgModel = ctauBkgModel = new RooAddPdf("pdfTot_Bkg","pdfTot_Bkg",*ws->pdf("pdfCTAUCOND_Bkg"), *ws->var("N_Bkg"));
	pdfTot_Bkg = new RooAddPdf("pdfTot_Bkg","",RooArgList(*ws->pdf("pdfCTAUCOND_Bkg")),RooArgList(*ws->var("N_Bkg")));
	//ws->import(*pdfCTAUCOND_Bkg);

	TCanvas* c_E =  new TCanvas("canvas_E","My plots",554,565,550,520);
	c_E->cd();
	TPad *pad_E_1 = new TPad("pad_E_1", "pad_E_1", 0, 0.16, 0.98, 1.0);
	pad_E_1->SetTicks(1,1);
	pad_E_1->Draw(); pad_E_1->cd();
	RooPlot* myPlot_E = ws->var("ctau3D")->frame(Bins(nCtauBins), Range(ctauLow, ctauHigh)); // bins
	myPlot_E->SetTitle("");

	pad_E_1->cd();
	gPad->SetLogy();
	RooPlot* myPlot2_E = (RooPlot*)myPlot_E->Clone();
	RooFitResult* fitCtauBkg = ws->pdf("pdfTot_Bkg")->fitTo(*dataw_Bkg, Save(), Range("ctauRange"), Extended(kTRUE), NumCPU(10), PrintLevel(-1));
	ws->data("dataw_Bkg")->plotOn(myPlot2_E,Name("data_ctauBkg"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kBlue+2), MarkerColor(kBlue+2));
	ws->pdf("pdfTot_Bkg")->plotOn(myPlot2_E,Name("pdfCTAUCOND_Bkg"),
			Normalization(ws->data("dataw_Bkg")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
			ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE),
			FillStyle(1001), FillColor(kAzure-9), LineColor(kBlack), LineWidth(3), DrawOption("LCF"), Precision(1e-4));
	ws->pdf("pdfTot_Bkg")->plotOn(myPlot2_E,Name("pdfCTAUCOND_BkgPR"),
			Normalization(ws->data("dataw_Bkg")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
			ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Bkg"), kTRUE), Components("pdfCTAUCOND_BkgPR"),
			LineWidth(2), LineColor(kRed+2));
	ws->data("dataw_Bkg")->plotOn(myPlot2_E,Name("data_ctauBkg"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kBlue+2), MarkerColor(kBlue+2));
	myPlot2_E->GetYaxis()->SetRangeUser(10e-2, 10e7);
	myPlot2_E->GetXaxis()->SetRangeUser(-4, 6);
	myPlot2_E->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
	myPlot2_E->SetFillStyle(4000);
	myPlot2_E->GetYaxis()->SetTitleOffset(1.43);
	myPlot2_E->GetXaxis()->SetLabelSize(0);
	myPlot2_E->GetXaxis()->SetTitleSize(0);
	myPlot2_E->Draw();
	TLegend* leg_E = new TLegend(text_x+0.25,text_y+0.03,text_x+0.38,text_y-0.17); leg_E->SetTextSize(text_size);
	leg_E->SetTextFont(43);
	leg_E->SetBorderSize(0);
	leg_E->AddEntry(myPlot2_E->findObject("data_ctauBkg"),"Data_Bkg","pe");
	leg_E->AddEntry(myPlot2_E->findObject("pdfCTAUCOND_Bkg"),"Total PDF","fl");
	//leg_E->AddEntry(myPlot2_E->findObject("test"),"?? PDF","l");
	leg_E->Draw("same");
	drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
	if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
	else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
	drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);

	drawText(Form("N_{Bkg} = %.f #pm %.f", ws->var("N_Bkg")->getVal(), ws->var("N_Bkg")->getError() ),text_x+0.5,text_y,text_color,text_size);
	drawText(Form("b_{Bkg} = %.4f #pm %.4f", ws->var("b_Bkg")->getVal(), ws->var("b_Bkg")->getError() ),text_x+0.5,text_y-y_diff*1,text_color,text_size);
	drawText(Form("fDFSS = %.4f #pm %.4f", ws->var("fDFSS")->getVal(), ws->var("fDFSS")->getError() ),text_x+0.5,text_y-y_diff*2,text_color,text_size);
	drawText(Form("fDLIV = %.4f #pm %.4f", ws->var("fDLIV")->getVal(), ws->var("fDLIV")->getError() ),text_x+0.5,text_y-y_diff*3,text_color,text_size);
	drawText(Form("#lambdaDDS_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDDS_Bkg")->getVal(), ws->var("lambdaDDS_Bkg")->getError() ),text_x+0.5,text_y-y_diff*4,text_color,text_size);
	drawText(Form("#lambdaDF_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDF_Bkg")->getVal(), ws->var("lambdaDF_Bkg")->getError() ),text_x+0.5,text_y-y_diff*5,text_color,text_size);
	drawText(Form("#lambdaDSS_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDSS_Bkg")->getVal(), ws->var("lambdaDSS_Bkg")->getError() ),text_x+0.5,text_y-y_diff*6,text_color,text_size);

	TPad *pad_E_2 = new TPad("pad_E_2", "pad_E_2", 0, 0.006, 0.98, 0.227);
	RooPlot* frameTMP_E = (RooPlot*)myPlot2_E->Clone("TMP");
	RooHist* hpull_E;
	pullDist(ws, pad_E_2, c_E, frameTMP_E, hpull_E, "data_ctauBkg", "pdfCTAUCOND_Bkg", "ctau3D", nCtauBins, ctauLow, ctauHigh, "#font[12]{l}_{J/#psi} (mm)");
	printChi2_test(ws, pad_E_2, frameTMP_E, hpull_E, fitCtauBkg, "ctau3D", "data_ctauBkg", "pdfCTAUCOND_Bkg", nCtauBins);
	pad_E_2->Update();

	c_E->Update();
	c_E->SaveAs(Form("../figs/2Dfit/Bkg_%s.pdf",kineLabel.Data()));
	RooArgSet* fitargs = new RooArgSet();
	fitargs->add(fitCtauBkg->floatParsFinal());
	RooDataSet *datasetCBkg = new RooDataSet("datasetCBkg","dataset with Ctau Background Fit result", *fitargs);
	RooWorkspace *wscbkg = new RooWorkspace("workspaceCBkg");
	wscbkg->import(*fitCtauBkg);
	//wscbkg->import(*pdfCTAUCOND_Bkg);
	
	ws->Print();

	TFile *outFile = new TFile(Form("../roots/2DFit/CtauBkgResult_%s.root",kineLabel.Data()),"recreate");
	fitCtauBkg->Write();
	pdfTot_Bkg->Write();
	//pdfCTAUCOND_Bkg->Write();
	//pdfTot_Bkg->Write();
	datasetCBkg->Write();
	wscbkg->Write();
	outFile->Close();


}
