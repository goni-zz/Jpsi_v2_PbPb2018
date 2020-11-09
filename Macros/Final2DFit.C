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




}



