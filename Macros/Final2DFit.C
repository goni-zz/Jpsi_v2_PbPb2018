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

void Final2DFit(
		float ptLow=3, float ptHigh=6.5,
		float yLow=1.6, float yHigh=2.4,
		int cLow=0, int cHigh=200,
		float muPtCut=0.0,
		bool whichModel=0,  // Nominal=0, Alternative=1
		int ICset=1
		)
{
	gStyle->SetEndErrorSize(0);
	TFile* f1; TFile* f2; TFile* f3; TFile* fMass; TFile* fCErr; TFile* fCRes; TFile* fCBkg; TFile* fCTrue;
	TString kineCut;
	TString SigCut;
	TString BkgCut;
	TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

	f1 = new TFile("../skimmedFiles/OniaRooDataSet_isMC0_JPsi1SW_2020819.root");
	f2 = new TFile("../skimmedFiles/OniaRooDataSet_JPsi_GENONLY_NonPrompt_20201006_365_2_ver618.root");
	f3 = new TFile("../skimmedFiles/OniaRooDataSet_JPsi_GENONLY_NonPrompt_20201006_6550_ver618.root");
	fMass = new TFile(Form("../roots/2DFit/MassFitResult_%s.root",kineLabel.Data()));
	fCErr = new TFile(Form("../roots/2DFit/CtauErrResult_%s.root",kineLabel.Data()));
	fCRes = new TFile(Form("../roots/2DFit/CtauResResult_%s.root",kineLabel.Data()));
	fCBkg = new TFile(Form("../roots/2DFit/CtauBkgResult_%s.root",kineLabel.Data()));
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
	RooAddPdf* pdfTot_Bkg = (RooAddPdf*)fCBkg->Get("pdfTot_Bkg");
	RooAddPdf* TrueModel_Tot = (RooAddPdf*)fCTrue->Get("TrueModel_Tot");
	RooDataSet *datasetMC = (RooDataSet*)f2->Get("dataset");
	RooHistPdf* pdfCTAUERR_Tot = (RooHistPdf*)fCErr->Get("pdfCTAUERR_Tot");
	RooHistPdf* pdfCTAUERR_Jpsi = (RooHistPdf*)fCErr->Get("pdfCTAUERR_Jpsi");
	RooHistPdf* pdfCTAUERR_Bkg = (RooHistPdf*)fCErr->Get("pdfCTAUERR_Bkg");
	if(ptLow >= 3.0 && ptHigh <= 6.5) {
		datasetMC = (RooDataSet*)f2->Get("dataset");
	}
	else {
		datasetMC = (RooDataSet*)f3->Get("dataset");}
	RooWorkspace *ws = new RooWorkspace("workspace");
    ws->import(*dataset);
    ws->import(*datasetMass);
    ws->import(*pdfMASS_Tot);
    ws->import(*dataw_Bkg);
    ws->import(*dataw_Sig);
	ws->import(*GaussModel_Tot);
	ws->import(*pdfTot_Bkg);
	ws->import(*TrueModel_Tot);
	ws->import(*pdfCTAUERR_Tot);
	ws->import(*pdfCTAUERR_Jpsi);
	ws->import(*pdfCTAUERR_Bkg);
	cout << "####################################" << endl;
	RooDataSet *reducedDS_A = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), SigCut.Data() );
	RooDataSet *reducedDS_B = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), BkgCut.Data() );
	reducedDS_A->SetName("reducedDS_A");
	reducedDS_B->SetName("reducedDS_B");
	//ws->import(*reducedDS_A);
	//ws->import(*reducedDS_B);
	//ws->var("mass")->setRange(massLow, massHigh);
	//ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
	//ws->var("ctau3DErr")->setRange(ctauLow, ctauHigh);
	//ws->var("ctau3DRes")->setRange(ctauResLow, ctauResHigh);
	//reducedDS_A->Print();
	//reducedDS_B->Print();
	RooCategory tp("tp","tp");
	tp.defineType("A");
	tp.defineType("B");

	// Create a dataset that imports contents of all the above datasets mapped by index category tp
	//RooDataSet* dsAB = new RooDataSet("dsAB","dsAB",RooArgSet(*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))),Index(tp),Import("A",*reducedDS_A),Import("B",*reducedDS_B));
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

	RooWorkspace *wsmc = new RooWorkspace("workspaceMC");
	wsmc->import(*datasetMC);
	tp.defineType("C");
	RooDataSet* reducedDS_MC = new RooDataSet("reducedDS_MC","reducedDS_MC",RooArgSet(*(wsmc->var("ctau3Dtrue")), *(wsmc->var("mass")), *(wsmc->var("pt")), *(wsmc->var("y"))),Index(tp),Import("C",*datasetMC));
	reducedDS_MC->SetName("reducedDS_MC");
	reducedDS_MC->Print();
	ws->import(*reducedDS_MC);
	//ws->var("ctau3Dtrue")->setRange(0, 5);
	ws->var("ctau3Dtrue")->setRange("ctauTrueRange", 0, 6.);
	ws->Print();


	//***********************************************************************
	//*************************** DRAW CTAU FIT *****************************
	//***********************************************************************
	ws->var("lambdaDSS")->setConstant(kTRUE);
	//make jpsi pdf
	ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUCOND_JpsiNoPR", "ctau3D",
				"lambdaDSS", "pdfCTAURES")); //NP
	ws->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_JpsiPR",//PR
				"pdfCTAURES"//resolution model
				));
	ws->factory("b_Jpsi[0.30, 0., 1.]");//NP fraction for Sig

	RooProdPdf pdfbkgPR("pdfCTAU_BkgPR", "", *ws->pdf("pdfCTAUERR_Bkg"),
			Conditional( *ws->pdf("pdfCTAUCOND_BkgPR"), RooArgList(*ws->var("ctau3D"))));
	ws->import(pdfbkgPR);
	RooProdPdf pdfbkgNoPR("pdfCTAU_BkgNoPR", "", *ws->pdf("pdfCTAUERR_Bkg"),
			Conditional( *ws->pdf("pdfCTAUCOND_BkgNoPR"), RooArgList(*ws->var("ctau3D"))));
	ws->import(pdfbkgNoPR);
	ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgPR",
				"pdfCTAU_BkgPR",
				"pdfMASS_bkg"
				));
	ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgNoPR",
				"pdfCTAU_BkgNoPR",
				"pdfMASS_bkg"
				));
	ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUMASS_Bkg",
				"b_Bkg",
				"pdfCTAUMASS_BkgNoPR",
				"pdfCTAUMASS_BkgPR"
				));
	RooProdPdf pdfJpsiPR("pdfCTAU_JpsiPR","",*ws->pdf("pdfCTAUERR_Jpsi"),
			Conditional(*ws->pdf("pdfCTAUCOND_JpsiPR"),RooArgList(*ws->var("ctau3D"))));
	ws->import(pdfJpsiPR);
	RooProdPdf pdfJpsiNoPR("pdfCTAU_JpsiNoPR","",*ws->pdf("pdfCTAUERR_Jpsi"),
			Conditional(*ws->pdf("pdfCTAUCOND_JpsiNoPR"),RooArgList(*ws->var("ctau3D"))));
	ws->import(pdfJpsiNoPR);

	ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_JpsiPR",
				"pdfCTAU_JpsiPR",
				"pdfMASS_Jpsi"
				));
	ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_JpsiNoPR",
				"pdfCTAU_JpsiNoPR",
				"pdfMASS_Jpsi"
				));
	ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUMASS_Jpsi",
				"b_Jpsi",
				"pdfCTAUMASS_JpsiNoPR",
				"pdfCTAUMASS_JpsiPR"
				));
	RooAbsPdf *themodel =NULL;
	themodel = new RooAddPdf("pdfCTAUMASS_Tot", "pdfCTAUMASS_Tot",
			RooArgList(*ws->pdf("pdfCTAUMASS_Bkg"), *ws->pdf("pdfCTAUMASS_Jpsi")),
			RooArgList(*ws->var("N_Bkg"), *ws->var("N_Jpsi")) );
	ws->import(*themodel);
	ws->pdf("pdfMASS_Tot")->getParameters(RooArgSet(*ws->var("mass")))->setAttribAll("Constant", kTRUE);
	std::vector< std::string > objs = {"Bkg", "Jpsi"};
	RooArgSet pdfList = RooArgSet("ConstraionPdfList");
	for (auto obj : objs) {
		if (ws->var(Form("N_%s", obj.c_str())))  {
			ws->factory(Form("Gaussian::%s_Gauss(%s,%s_Mean[%f],%s_Sigma[%f])",
						Form("N_%s", obj.c_str()), Form("N_%s", obj.c_str()),
						Form("N_%s", obj.c_str()), ws->var(Form("N_%s", obj.c_str()))->getValV(),
						Form("N_%s", obj.c_str()), ws->var(Form("N_%s", obj.c_str()))->getError()));

			pdfList.add(*ws->pdf(Form("N_%s_Gauss", obj.c_str())), kFALSE);
			std::cout << "[INFO] Constraining N_" << obj << " with Mean : " << ws->var(Form("N_%s_Mean", obj.c_str()))->getVal()
				<< " and Sigma: " << ws->var(Form("N_%s_Sigma", obj.c_str()))->getVal() << std::endl;
		}
	}
	ws->defineSet("ConstrainPdfList", pdfList);
	ws->pdf("pdfCTAURES")->getParameters(RooArgSet(*ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("ctau3DErr")))->setAttribAll("Constant", kTRUE);
	ws->pdf("pdfCTAU_JpsiPR")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr")))->setAttribAll("Constant", kTRUE);
	ws->pdf("pdfCTAU_BkgPR")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr")))->setAttribAll("Constant", kTRUE);
	ws->pdf("pdfCTAU_BkgNoPR")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr")))->setAttribAll("Constant", kTRUE);

	RooArgSet* params = (RooArgSet*) ws->pdf("pdfCTAUMASS_Tot")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("mass"), *ws->var("ctau3DErr")));
	ws->saveSnapshot(("pdfCTAUMASS_Tot_parIni"),*params,kTRUE);
	delete params;

	RooArgSet *newpars = (RooArgSet*) ws->pdf("pdfCTAUMASS_Tot")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr"), *ws->var("ctau3DRes"), *ws->var("mass")));

	TCanvas* c_G =  new TCanvas("canvas_G","My plots",1108,565,550,520);
	c_G->cd();
	TPad *pad_G_1 = new TPad("pad_G_1", "pad_G_1", 0, 0.16, 0.98, 1.0);
	pad_G_1->SetTicks(1,1);
	pad_G_1->Draw(); pad_G_1->cd();
	RooPlot* myPlot_G = ws->var("ctau3D")->frame(nCtauBins); // bins
	myPlot_G->SetTitle("");

	c_G->cd();
	c_G->SetLogy();
	double normDSTot = ws->data("dsAB")->sumEntries()/ws->data("dsAB")->sumEntries();
	double normBkg = ws->data("dsAB")->sumEntries()*normDSTot/ws->data("dataw_Bkg")->sumEntries();
	double normJpsi =ws->data("dsAB")->sumEntries()*normDSTot/ws->data("dataw_Sig")->sumEntries();

	cout<<"##############START TOTAL CTAU FIT############"<<endl;
	RooFitResult* fitResult = ws->pdf("pdfCTAUMASS_Tot")->fitTo(*dsAB, Extended(kTRUE), ExternalConstraints(*ws->set("ConstrainPdfList")),
			NumCPU(12), SumW2Error(true), PrintLevel(-1), Save());
	ws->import(*fitResult, "fitResult_pdfCTAUMASS_Tot");

	//DRAW
	RooPlot* myPlot2_G = (RooPlot*)myPlot_G->Clone();
	myPlot2_G->updateNormVars(RooArgSet(*ws->var("mass"), *ws->var("ctau3D"), *ws->var("ctau3DErr")));
	ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G,Name("PDF"),
			ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsAB"), kTRUE),
			Normalization(normDSTot, RooAbsReal::NumEvent),
			FillStyle(1001), FillColor(kViolet+6), VLines(), DrawOption("LF"), NumCPU(32), LineColor(kBlack)
			);
	ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G,Name("BKG"),Components(RooArgSet( *ws->pdf("pdfCTAUMASS_Bkg") )),
			ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsAB"), kTRUE),
			Normalization(normDSTot, RooAbsReal::NumEvent),
			FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LF"), NumCPU(32)
			);
	ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G,Name("JPSIPR"),Components(RooArgSet( *ws->pdf("pdfCTAUMASS_JpsiPR") )),
			ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsAB"), kTRUE),
			Normalization(normDSTot, RooAbsReal::NumEvent),
			LineColor(kRed+3), Precision(1e-5), NumCPU(32)
			);
	ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G,Name("JPSINOPR"),Components(RooArgSet( *ws->pdf("pdfCTAUMASS_JpsiNoPR"))),
			ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsAB"), kTRUE),
			Normalization(normDSTot, RooAbsReal::NumEvent),
			LineColor(kGreen+3), Precision(1e-5), NumCPU(32)
			);
	ws->data("dsAB")->plotOn(myPlot2_G,Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack));
	//ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G,Name("PDFLINE"),
	//                                     ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsAB"), kTRUE),
	//                                     Normalization(normDSTot, RooAbsReal::NumEvent),
	//                                     LineColor(kBlack), NumCPU(32)
	//                                     );
	ws->saveSnapshot("pdfCTAUMASS_Tot_parFit",*newpars,kTRUE);
	myPlot2_G->GetYaxis()->SetRangeUser(10e-2, 10e7);
	myPlot2_G->GetXaxis()->SetRangeUser(-4, 6);
	myPlot2_G->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
	myPlot2_G->Draw();
	TLegend* leg_G = new TLegend(text_x+0.25,text_y+0.03,text_x+0.38,text_y-0.17); leg_G->SetTextSize(text_size);
	leg_G->SetTextFont(43);
	leg_G->SetBorderSize(0);
	leg_G->AddEntry(myPlot2_G->findObject("dOS"),"Data","pe");
	leg_G->AddEntry(myPlot2_G->findObject("PDF"),"Total fit","fl");
	leg_G->AddEntry(myPlot2_G->findObject("BKG"),"Background","pe");
	leg_G->AddEntry(myPlot2_G->findObject("JPSIPR"),"J/#psi Prompt","l");
	leg_G->AddEntry(myPlot2_G->findObject("JPSINOPR"),"J/#psi Non-Prompt","l");
	leg_G->Draw("same");
	drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
	if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
	else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
	drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);

	c_G->Update();
	c_G->SaveAs(Form("../figs/2DFit/2DFit_%s.pdf",kineLabel.Data()));
	TString outFileName;
	if (whichModel){
		outFileName = Form("../roots/2DFit/test_Sim_BSplit_altfitresults_Jpsi_%s.root",kineLabel.Data());
	}
	else {
		outFileName = Form("../roots/2DFit/test_Sim_BSplit_nomfitresults_Jpsi_%s.root",kineLabel.Data());
	}
	TFile* outf = new TFile(outFileName,"recreate");
	ws->Write();
	outf->Close();

}



