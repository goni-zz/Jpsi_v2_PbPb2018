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
void CtauErr(
		float ptLow=3, float ptHigh=6.5,
		float yLow=1.6, float yHigh=2.4,
		int cLow=0, int cHigh=200,
		float muPtCut=0.0,
		bool whichModel=0,  // Nominal=0, Alternative=1
		int ICset=1
		)
{
	gStyle->SetEndErrorSize(0);
	TFile* f1; TFile* f2; TFile* f3; TFile* fMass;
	TString kineCut;
	TString SigCut;
	TString BkgCut;
	TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh);

	f1 = new TFile("../skimmedFiles/OniaRooDataSet_isMC0_JPsi1SW_2020819.root");
	fMass = new TFile(Form("../roots/2DFit/MassFitResult_%s.root",kineLabel.Data()));
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

	RooWorkspace *ws = new RooWorkspace("workspace");
	ws->import(*dataset);
	ws->import(*datasetMass);
	ws->import(*pdfMASS_Tot);
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
	//**************************** CTAU ERR FIT *****************************
	//***********************************************************************

	cout << endl << "*************** Start SPLOT *****************" << endl << endl;
	//SPlot Ctau Error
	RooRealVar *sigYield = ws->var("N_Jpsi");
	RooRealVar *bkgYield = ws->var("N_Bkg");
	RooArgList yieldList;
	yieldList.add(*ws->var("N_Jpsi"));
	yieldList.add(*ws->var("N_Bkg"));
	cout<<"Sig Yield: "<<sigYield->getVal()<<" +/- "<<sigYield->getError()<<endl;
	cout<<"Bkg Yield: "<<bkgYield->getVal()<<" +/- "<<bkgYield->getError()<<endl;
	RooDataSet* data = (RooDataSet*)ws->data("dsAB")->Clone("TMP_DATA");
	RooArgSet* cloneSet = (RooArgSet*)RooArgSet(*ws->pdf("pdfMASS_Tot"),"pdfMASS_Tot").snapshot(kTRUE);
	auto clone_mass_pdf = (RooAbsPdf*)cloneSet->find("pdfMASS_Tot");
	clone_mass_pdf->setOperMode(RooAbsArg::ADirty, kTRUE);

	RooStats::SPlot sData = RooStats::SPlot("sData", "An SPlot", *data, clone_mass_pdf, yieldList);
	ws->import(*data, Rename("dataset_SPLOT"));
	cout<<"[INFO] Jpsi yield -> Mass Fit:"<<ws->var("N_Jpsi")->getVal()<<", sWeights :"<<sData.GetYieldFromSWeight("N_Jpsi")<<endl;
	cout<<"[INFO] Bkg  yield -> Mass Fit:"<<ws->var("N_Bkg")->getVal()<<", sWeights :"<<sData.GetYieldFromSWeight("N_Bkg")<<endl;
	//create weighted data sets
	//total
	TH1D* hTot = (TH1D*)ws->data("dsAB")->createHistogram(("hTot"), *ws->var("ctau3DErr"),Binning(nCtauErrBins,ctauErrLow,ctauErrHigh));
	RooDataHist* totHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hTot);
	RooHistPdf* pdfCTAUERR_Tot = new RooHistPdf("pdfCTAUERR_Tot","hist pdf", *ws->var("ctau3DErr"), *totHist);
	//bkg
	RooDataSet* dataw_Bkg = new RooDataSet("dataw_Bkg","TMP_BKG_DATA", (RooDataSet*)ws->data("dataset_SPLOT"),
			RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Bkg_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D")), 0, "N_Bkg_sw");
	TH1D* hBkg = (TH1D*)dataw_Bkg->createHistogram(("hBkg"), *ws->var("ctau3DErr"),Binning(nCtauErrBins,ctauErrLow,ctauErrHigh));
	RooDataHist* bkgHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hBkg);
	RooHistPdf* pdfCTAUERR_Bkg = new RooHistPdf("pdfCTAUERR_Bkg","hist pdf", *ws->var("ctau3DErr"), *bkgHist);
	//data
	RooDataSet* dataw_Sig = new RooDataSet("dataw_Sig","TMP_SIG_DATA", (RooDataSet*)ws->data("dataset_SPLOT"),
			RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Jpsi_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D")), 0, "N_Jpsi_sw");
	TH1D* hSig = (TH1D*)dataw_Sig->createHistogram(("hSig"), *ws->var("ctau3DErr"),Binning(nCtauErrBins,ctauErrLow,ctauErrHigh));
	RooDataHist* sigHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hSig);
	RooHistPdf* pdfCTAUERR_Jpsi = new RooHistPdf("pdfCTAUERR_Jpsi","hist pdf", *ws->var("ctau3DErr"), *sigHist);
	//import
	ws->import(*dataw_Sig);
	ws->import(*dataw_Bkg);
	ws->import(*pdfCTAUERR_Tot);
	ws->import(*pdfCTAUERR_Jpsi);
	ws->import(*pdfCTAUERR_Bkg);
	cout<<ws->data("dsAB")->numEntries()<<endl;
	cout<<ws->data("dataset_SPLOT")->numEntries()<<endl;



	TCanvas* c_B =  new TCanvas("canvas_B","My plots",554,4,550,520);
	c_B->cd();
	TPad *pad_B_1 = new TPad("pad_B_1", "pad_B_1", 0, 0.16, 0.98, 1.0);
	pad_B_1->SetTicks(1,1);
	pad_B_1->Draw(); pad_B_1->cd();
	RooPlot* myPlot_B = ws->var("ctau3DErr")->frame(nCtauErrBins); // bins
	myPlot_B->SetTitle("");
	c_B->cd();
	c_B->SetLogy();

	pad_B_1->cd();
	gPad->SetLogy();
	RooPlot* myPlot2_B = (RooPlot*)myPlot_B->Clone();
	ws->data("dataset_SPLOT")->plotOn(myPlot2_B,Name("dataCTAUERR_Tot"), MarkerSize(.7), Binning(nCtauErrBins));
	ws->pdf("pdfCTAUERR_Tot")->plotOn(myPlot2_B,Name("pdfCTAUERR_Tot"), LineColor(kGreen+1), Range("ctauErrRange"), LineWidth(2));
	ws->data("dataw_Sig")->plotOn(myPlot2_B,Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2), Binning(nCtauErrBins));
	ws->pdf("pdfCTAUERR_Jpsi")->plotOn(myPlot2_B,Name("pdfCTAUERR_Jpsi"),LineColor(kRed+2), LineWidth(2), Range("ctauErrRange"));
	ws->data("dataw_Bkg")->plotOn(myPlot2_B,Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue+2), MarkerColor(kBlue+2), Binning(nCtauErrBins));
	ws->pdf("pdfCTAUERR_Bkg")->plotOn(myPlot2_B,Name("pdfCTAUERR_Bkg"), LineColor(kBlue+2), LineWidth(2), Range("ctauErrRange"));
	myPlot2_B->GetYaxis()->SetRangeUser(10e-4, 10e8);
	myPlot2_B->GetXaxis()->CenterTitle();
	myPlot2_B->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
	myPlot2_B->SetFillStyle(4000);
	myPlot2_B->GetYaxis()->SetTitleOffset(1.43);
	myPlot2_B->GetXaxis()->SetLabelSize(0);
	myPlot2_B->GetXaxis()->SetTitleSize(0);
	myPlot2_B->Draw();
	//Double_t outTot = ws->data("dsAB")->numEntries();
	//Double_t outErr = ws->data("dsAB")->reduce(Form("(ctauErr>=%.6f || ctauErr<=%.6f)", ctauErrHigh, ctauErrLow))->numEntries();
	//cout<<(outErr*100)/outTot<<endl;
	TLegend* leg_B = new TLegend(text_x+0.5,text_y-0.2,text_x+0.7,text_y); leg_B->SetTextSize(text_size);
	leg_B->SetTextFont(43);
	leg_B->SetBorderSize(0);
	leg_B->AddEntry(myPlot2_B->findObject("dataCTAUERR_Tot"),"Data","pe");
	leg_B->AddEntry(myPlot2_B->findObject("pdfCTAUERR_Tot"),"Total PDF","l");
	leg_B->AddEntry(myPlot2_B->findObject("pdfCTAUERR_Jpsi"),"Signal","l");
	leg_B->AddEntry(myPlot2_B->findObject("pdfCTAUERR_Bkg"),"Background","l");
	leg_B->Draw("same");
	drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",ptLow, ptHigh ),text_x,text_y,text_color,text_size);
	if(yLow==0)drawText(Form("|y^{#mu#mu}| < %.1f",yHigh), text_x,text_y-y_diff,text_color,text_size);
	else if(yLow!=0)drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh), text_x,text_y-y_diff,text_color,text_size);
	//drawText(Form("Cent. %d - %d%s", cLow, cHigh, "%"),text_x,text_y-y_diff*2,text_color,text_size);
	//drawText(Form("n_{J/#psi} = %.f #pm %.f",sData.GetYieldFromSWeight("N_Jpsi"),ws->var("N_Jpsi")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
	//drawText(Form("n_{Bkg} = %.f #pm %.f",sData.GetYieldFromSWeight("N_Bkg"), ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
	TPad *pad_B_2 = new TPad("pad_B_2", "pad_B_2", 0, 0.006, 0.98, 0.227);
	c_B->cd();
	pad_B_2->Draw();
	pad_B_2->cd();
	pad_B_2->SetTopMargin(0); // Upper and lower plot are joined
	pad_B_2->SetBottomMargin(0.67);
	pad_B_2->SetBottomMargin(0.4);
	pad_B_2->SetFillStyle(4000);
	pad_B_2->SetFrameFillStyle(4000);
	pad_B_2->SetTicks(1,1);

	RooPlot* frameTMP_B = (RooPlot*)myPlot2_B->Clone("TMP");
	RooHist* hpull_B = frameTMP_B->pullHist("dataCTAUERR_Tot","pdfCTAUERR_Tot");
	hpull_B->SetMarkerSize(0.8);
	RooPlot* pullFrame_B = ws->var("ctau3DErr")->frame(Title("Pull Distribution")) ;
	pullFrame_B->addPlotable(hpull_B,"PX") ;
	pullFrame_B->SetTitle("");
	pullFrame_B->SetTitleSize(0);
	pullFrame_B->GetYaxis()->SetTitleOffset(0.3) ;
	pullFrame_B->GetYaxis()->SetTitle("Pull") ;
	pullFrame_B->GetYaxis()->SetTitleSize(0.15) ;
	pullFrame_B->GetYaxis()->SetLabelSize(0.15) ;
	pullFrame_B->GetYaxis()->SetRangeUser(-3.8,3.8);
	pullFrame_B->GetYaxis()->CenterTitle();

	pullFrame_B->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
	pullFrame_B->GetXaxis()->SetTitleOffset(1.05) ;
	pullFrame_B->GetXaxis()->SetLabelOffset(0.04) ;
	pullFrame_B->GetXaxis()->SetLabelSize(0.15) ;
	pullFrame_B->GetXaxis()->SetTitleSize(0.15) ;
	pullFrame_B->GetXaxis()->CenterTitle();

	pullFrame_B->GetYaxis()->SetTickSize(0.04);
	pullFrame_B->GetYaxis()->SetNdivisions(404);
	pullFrame_B->GetXaxis()->SetTickSize(0.03);
	pullFrame_B->Draw() ;

	TLine *lB = new TLine(ctauErrLow,0, ctauErrHigh,0);
	lB->SetLineStyle(1);
	lB->Draw("same");
	pad_B_2->Update();
	cout << endl << "************** Finished SPLOT *****************" << endl << endl;

	c_B->Update();
	c_B->SaveAs(Form("../figs/2DFit/ctauErr_%s.pdf",kineLabel.Data()));
	c_B->SaveAs(Form("../figs/2DFit/ctauErr_%s.png",kineLabel.Data()));


	TFile *outFile = new TFile(Form("../roots/2DFit/CtauErrResult_%s.root",kineLabel.Data()),"recreate");
	dataw_Bkg->Write();
	dataw_Sig->Write();
	pdfCTAUERR_Tot->Write();
	pdfCTAUERR_Bkg->Write();
	pdfCTAUERR_Jpsi->Write();
	outFile->Close();


}
