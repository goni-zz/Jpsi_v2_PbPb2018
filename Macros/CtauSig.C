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

void CtauSig(
		float ptLow=3, float ptHigh=6.5,
		float yLow=1.6, float yHigh=2.4,
		int cLow=0, int cHigh=200,
		float muPtCut=0.0,
		bool whichModel=0,  // Nominal=0, Alternative=1
		int ICset=1
		)
{
	gStyle->SetEndErrorSize(0);
	TFile* f1; TFile* fMass; TFile* fCErr; TFile* fCRes; TFile* fCTrue;
	TString kineCut;
	TString SigCut;
	TString BkgCut;
	TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

	f1 = new TFile("../skimmedFiles/OniaRooDataSet_isMC0_JPsi1SW_2020819.root");
	fMass = new TFile(Form("../roots/2DFit/MassFitResult_%s.root",kineLabel.Data()));
	fCErr = new TFile(Form("../roots/2DFit/CtauErrResult_%s.root",kineLabel.Data()));
	fCRes = new TFile(Form("../roots/2DFit/CtauResResult_%s.root",kineLabel.Data()));
	fCTrue = new TFile(Form("../roots/2DFit/CtauTrueResult_%s.root",kineLabel.Data()));

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
	RooAddPdf* TrueModel_Tot = (RooAddPdf*)fCTrue->Get("TrueModel_Tot");

	RooWorkspace *ws = new RooWorkspace("workspace");
	ws->import(*dataset);
	ws->import(*datasetMass);
	ws->import(*pdfMASS_Tot);
	ws->import(*dataw_Bkg);
	ws->import(*dataw_Sig);
	ws->import(*GaussModel_Tot);
	ws->import(*TrueModel_Tot);
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
	//**************************** Sig CTAU FIT *****************************
	//***********************************************************************
	cout << endl << "************** Start SIG Ctau Fit *****************" << endl << endl;
	//make parameter 3 exp
	ws->factory("zeroMean[0.]");
	ws->factory("fDFSS[0.5, 1e-6, 1.1]");
	ws->factory("fDLIV[0.5, 1e-6, 1.]");
	ws->factory("fDFSS_Sig[0.5, 1e-6, 1.]");
	ws->factory("fDLIV_Sig[0.5, 1e-6, 1.]");
	ws->factory("lambdaDDS_Sig[0.5, 1e-6, 1.]");
	ws->factory("lambdaDF_Sig[0.5, 1e-6, 1.]");
	//make parameter exp-->should change as an initial
	//parameters fixed by Resolution model
	ws->var("ctauRes_mean")->setConstant(kTRUE);
	ws->var("s1_CtauRes")->setConstant(kTRUE);
	ws->var("s2_CtauRes")->setConstant(kTRUE);
	ws->var("s3_CtauRes")->setConstant(kTRUE);
	ws->var("f2_CtauRes")->setConstant(kTRUE);
	ws->var("f_CtauRes")->setConstant(kTRUE);
	ws->var("zeroMean")->setConstant(kTRUE);

	ws->var("lambdaDSS")->setConstant(kTRUE);
	//make res model
	ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauResSig1", "ctau3D",
				"ctauRes_mean", //"ctau1_CtauRes",
				"s1_CtauRes",
				"zeroMean",
				"ctau3DErr"
				));
	ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauResSig2", "ctau3D",
				"ctauRes_mean", //"ctau2_CtauRes",
				"s2_CtauRes",
				"zeroMean",
				"ctau3DErr"
				));
	ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauResSig3", "ctau3D",
				"ctauRes_mean", //"ctau3_CtauRes",
				"s3_CtauRes",
				"zeroMean",
				"ctau3DErr"
				));
	ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "ctauResSig32", "ctauResSig3", "ctauResSig2", "f2_CtauRes"));
	ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "pdfSigCTAURES", "ctauResSig32", "ctauResSig1", "f_CtauRes"));
	//make 3 exp
	ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUCOND_SigNoPR", "ctau3D",
				"lambdaDSS", "pdfSigCTAURES"));

	ws->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_SigPR",//PR
				"pdfSigCTAURES"//resolution model
				));
	ws->factory("b_Frac[0.30, 0., 1.]");//NP fraction for Sig
	ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_Sig",
				"b_Frac",
				"pdfCTAUCOND_SigNoPR",
				"pdfCTAUCOND_SigPR"
				));
	ws->factory(Form("RooExtendPdf::%s(%s,%s)", "pdfTot_Sig", "pdfCTAUCOND_Sig", "N_Jpsi")); //N_Jpsi is number of bkg from dataw_Sig
	RooAddPdf* pdfTot_Sig = new RooAddPdf("pdfTot_Sig","",RooArgList(*ws->pdf("pdfCTAUCOND_Sig")), RooArgList(*ws->pdf("N_Sig")));

	TCanvas* c_F =  new TCanvas("canvas_F","My plots",1108,565,550,520);
	c_F->cd();
	TPad *pad_F_1 = new TPad("pad_F_1", "pad_F_1", 0, 0.16, 0.98, 1.0);
	pad_F_1->SetTicks(1,1);
	pad_F_1->Draw(); pad_F_1->cd();
	RooPlot* myPlot_F = ws->var("ctau3D")->frame(nCtauBins); // bins
	myPlot_F->SetTitle("");

	pad_F_1->cd();
	gPad->SetLogy();
	RooPlot* myPlot2_F = (RooPlot*)myPlot_F->Clone();
	cout << "HERE" << endl;
	RooFitResult* fitCtauSig = ws->pdf("pdfTot_Sig")->fitTo(*dataw_Sig, Save(), Range("ctauRange"), Extended(kTRUE), NumCPU(12), PrintLevel(-1));
	cout << "HERE!" << endl;
	ws->data("dataw_Sig")->plotOn(myPlot2_F,Name("dataHist_ctauSig"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kRed+2), MarkerColor(kRed+2));
	ws->pdf("pdfTot_Sig")->plotOn(myPlot2_F, Name("ctauSig_Tot"),
			Normalization(ws->data("dataw_Sig")->sumEntries(), RooAbsReal::NumEvent), NormRange("catuRange"),
			ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Sig"), kTRUE),
			FillStyle(1001), FillColor(kRed-10), LineColor(kBlack), LineWidth(3), DrawOption("LCF"), Precision(1e-4));
	ws->pdf("pdfTot_Sig")->plotOn(myPlot2_F,Name("PR"),
			Normalization(ws->data("dataw_Sig")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
			ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Sig"), kTRUE), Components("pdfCTAUCOND_SigPR"),
			LineWidth(2), LineColor(kMagenta+3));
	ws->pdf("pdfTot_Sig")->plotOn(myPlot2_F,Name("NP"),
			Normalization(ws->data("dataw_Sig")->sumEntries(), RooAbsReal::NumEvent), NormRange("ctauRange"),
			ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataw_Sig"), kTRUE), Components("pdfCTAUCOND_SigNoPR"),
			LineWidth(2), LineColor(kOrange-5));
	ws->data("dataw_Sig")->plotOn(myPlot2_F,Name("dataHist_ctauSig"), DataError(RooAbsData::SumW2), MarkerSize(.7), Binning(nCtauBins), LineColor(kRed+2), MarkerColor(kRed+2));
	myPlot2_F->GetYaxis()->SetRangeUser(10e-2, 10e6);
	myPlot2_F->GetXaxis()->SetRangeUser(-4, 6);
	myPlot2_F->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
	myPlot2_F->SetFillStyle(4000);
	myPlot2_F->GetYaxis()->SetTitleOffset(1.43);
	myPlot2_F->GetXaxis()->SetLabelSize(0);
	myPlot2_F->GetXaxis()->SetTitleSize(0);
	myPlot2_F->Draw();

	TLegend* leg_F = new TLegend(text_x+0.25,text_y+0.03,text_x+0.37,text_y-0.17); leg_F->SetTextSize(text_size);
	leg_F->SetTextFont(43);
	leg_F->SetBorderSize(0);
	leg_F->AddEntry(myPlot2_F->findObject("dataHist_ctauSig"),"Data_Sig","pe");
	leg_F->AddEntry(myPlot2_F->findObject("ctauSig_Tot"),"Total fit","fl");
	leg_F->AddEntry(myPlot2_F->findObject("PR"),"Prompt J/#psi","l");
	leg_F->AddEntry(myPlot2_F->findObject("NP"),"Non-Prompt J/#psi","l");
	//leg_F->AddEntry(myPlot_F->findObject("test"),"?? PDF","l");
	leg_F->Draw("same");
	drawText(Form("%.1f < p_{T}^{#mu#mu} < %1.f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
	if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
	else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
	drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);

	drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError() ),text_x+0.5,text_y,text_color,text_size);
	drawText(Form("b_{J/#psi} = %.4f #pm %.4f", ws->var("b_Frac")->getVal(),ws->var("b_Frac")->getError() ),text_x+0.5,text_y-y_diff*1,text_color,text_size);
	drawText(Form("#lambdaDSS_{J/#psi} = %.4f #pm %.4f", ws->var("lambdaDSS")->getVal(), ws->var("lambdaDSS")->getError() ),text_x+0.5,text_y-y_diff*2,text_color,text_size);

	TPad *pad_F_2 = new TPad("pad_F_2", "pad_F_2", 0, 0.006, 0.98, 0.227);
	RooPlot* frameTMP_F = (RooPlot*)myPlot2_F->Clone("TMP");
	RooHist* hpull_F;
	pullDist(ws, pad_F_2, c_F, frameTMP_F, hpull_F, "dataHist_ctauSig", "ctauSig_Tot", "ctau3D", nCtauBins, ctauLow, ctauHigh, "#font[12]{l}_{J/#psi} (mm)");
	printChi2_test(ws, pad_F_2, frameTMP_F, hpull_F, fitCtauSig, "ctau3D", "dataHist_ctauSig", "ctauSig_Tot", nCtauBins);
	pad_F_2->Update();

	c_F->Update();
	c_F->SaveAs(Form("../figs/2DFit/Sig_%s.pdf",kineLabel.Data()));

	ws->Print();
	
	TFile *outFile = new TFile(Form("../roots/2DFit/CtauSigResult_%s.root",kineLabel.Data()),"recreate");
	pdfTot_Sig->Write();
	outFile->Close();

}
