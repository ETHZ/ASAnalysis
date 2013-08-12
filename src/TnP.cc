#include "TnP.hh"

#include "helper/AnaClass.hh"
#include "helper/Utilities.hh"
#include "helper/Monitor.hh"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TRandom3.h"

#include "TWbox.h"
#include "TMultiGraph.h"
#include "TGaxis.h"


#include "TF1.h"
#include "TH1.h"
#include "TList.h"
#include "TCanvas.h"
#include "TStyle.h"

// stuff for fitting
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
//#include "RooDoubleCB.h"
#include "RooBreitWigner.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"

#include "Linkdef.h"


#include <iostream>
#include <iomanip>
#include <time.h> // access to date/time

#include <typeinfo>


using namespace std;
using namespace RooFit;

//____________________________________________________________________________
TnP::TnP(TString inputfile, bool createHistos){
	fCreateHistos = createHistos;
	fInputFile = inputfile;

	fOutputFileName = "tnp_output.root";

	checkFlavor();

	if(fIsMu){
		fIsoCut = 0.10;
		fD0Cut  = 0.005;
	}
	else{
		fIsoCut = 0.09;
		fD0Cut  = 0.01;
	}
	
	effs = new eff[fnBins];
	// dynamic creation of the arrays
	fPassIDoverAll   = new TH1F*[fnBins];
	fPassIPoverID    = new TH1F*[fnBins];
	fPassISOoverIDIP = new TH1F*[fnBins];

	fFailIDoverAll   = new TH1F*[fnBins];
	fFailIPoverID    = new TH1F*[fnBins];
	fFailISOoverIDIP = new TH1F*[fnBins];

	if(fCreateHistos){
		for(int i=0; i<fnBins; ++i){
			fPassIDoverAll   [i] = new TH1F(Form("PassIDoverAll_%d"  ,i), Form("PassIDoverAll_%d"  ,i), 120, 0., 120.);
			fPassIPoverID    [i] = new TH1F(Form("PassIPoverID_%d"   ,i), Form("PassIPoverID_%d"   ,i), 120, 0., 120.);
			fPassISOoverIDIP [i] = new TH1F(Form("PassISOoverIDIP_%d",i), Form("PassISOoverIDIP_%d",i), 120, 0., 120.);
                                                                                                  
			fFailIDoverAll   [i] = new TH1F(Form("FailIDoverAll_%d"  ,i), Form("FailIDoverAll_%d"  ,i), 120, 0., 120.);
			fFailIPoverID    [i] = new TH1F(Form("FailIPoverID_%d"   ,i), Form("FailIPoverID_%d"   ,i), 120, 0., 120.);
			fFailISOoverIDIP [i] = new TH1F(Form("FailISOoverIDIP_%d",i), Form("FailISOoverIDIP_%d",i), 120, 0., 120.);

			fPassIDoverAll   [i] ->Sumw2();
			fPassIPoverID    [i] ->Sumw2();
			fPassISOoverIDIP [i] ->Sumw2();

			fFailIDoverAll   [i] ->Sumw2();
			fFailIPoverID    [i] ->Sumw2();
			fFailISOoverIDIP [i] ->Sumw2();
		}
	} // end if createhistos
	else{
		TFile *inFile = TFile::Open(fInputFile, "READ");
		if(inFile == NULL){
			cout << "File " << fInputFile << " does not exist!" << endl;
			exit(1);
		}
		else{
			std::cout << "reading file " + fInputFile << std::endl;
		}
		inFile->cd();

		for(int i=0; i<fnBins; ++i){
			if(fVerbose > 2) inFile->ls();
			getObjectSafe(inFile, Form("PassIDoverAll_%d"  ,i), fPassIDoverAll   [i]);
			getObjectSafe(inFile, Form("PassIPoverID_%d"   ,i), fPassIPoverID    [i]);
			getObjectSafe(inFile, Form("PassISOoverIDIP_%d",i), fPassISOoverIDIP [i]);

			getObjectSafe(inFile, Form("FailIDoverAll_%d"  ,i), fFailIDoverAll   [i]);
			getObjectSafe(inFile, Form("FailIPoverID_%d"   ,i), fFailIPoverID    [i]);
			getObjectSafe(inFile, Form("FailISOoverIDIP_%d",i), fFailISOoverIDIP [i]);
		}
	} // end else createhistos
}

void TnP::checkFlavor(){
	// check in the first event of the tree whether it's muon or electron tree
	TFile *pFile = TFile::Open(fInputFile);
	TTree *tree; getObjectSafe(pFile, "probeTree", tree);
	int isMu;
	tree->SetBranchAddress("isMuEvent", &isMu);
	tree->GetEntry(1);
	if(isMu == 1) {
		fIsMu  = true;
		fnBins = 13; // set the proper nbins
	}
	else {
		fIsMu  = false;
		fnBins = 25; // set the proper nbins
	}
	TString flav = fIsMu ? "muons" : "electrons";
	cout << "We're dealing with "+flav+" here!!" << endl;

	pFile->Close();
	delete pFile, tree;

}

//____________________________________________________________________________
TnP::~TnP(){
}

void TnP::doFitting(){

	if(fCreateHistos) {
		cout << "i'm going to do create all the histograms first, that will take a while... " << endl;
		fillHistos();
		writeHistos();
	}

	for(int i=0; i<fnBins; ++i){
		if(i>0) continue;
		simFitPassFail(fPassIDoverAll  [i], fFailIDoverAll  [i], 0, i); // ID
		simFitPassFail(fPassIPoverID   [i], fFailIPoverID   [i], 1, i); // IP
		simFitPassFail(fPassISOoverIDIP[i], fFailISOoverIDIP[i], 2, i); // ISO
	}

	// if (fIsMu) printMuTable();
	// else       printElTable();

}
void TnP::simFitPassFail(TH1F* passHisto, TH1F* failHisto, int flag, int bin){

	if(passHisto == NULL) cout << "histogram doesn't exist!!!!!" << endl;

	// passing probes
	RooRealVar px("px", "px", 0, 120) ;
	RooRealVar pmean ("pmean" , "pass mean of gaussian" , 90, 85 , 95);
	RooRealVar psigma("psigma", "pass width of gaussian",  1, 0. ,  5);

	RooRealVar a("a","a",3.,2.,10.);
    RooRealVar aDx("aDx","aDx",3.,2.,10.);
    RooRealVar n("n","n",5.,0.,10.);   
    RooRealVar nDx("nDx","nDx",5.,0.,10.);   

    RooDoubleCB func1("cb","cb PDF", px, pmean, psigma, a, n, aDx, nDx) ;
	// RooGaussian fgauss("pgauss", "pgaussian PDF", px, pmean, psigma) ;

	RooDataHist passData("passData", "passData", px, Import(*passHisto) );
	RooPlot* passFrame = px.frame(Title("pass TH1 with Poisson error bars"));
	passData.plotOn(passFrame);
	RooGaussian pgauss("pgauss", "pgaussian PDF", px, pmean, psigma) ;

	func1.fitTo (passData );
	func1.plotOn(passFrame);


//	// failing probes
//	RooRealVar fx("fx", "fx", 0, 120) ;
//	RooRealVar fmean ("fmean" , "fail mean of gaussian" , 90, 85 , 95);
//	RooRealVar fsigma("fsigma", "fail width of gaussian",  1, 0. ,  5);
//
//	RooDataHist failData("failData", "failData", px, Import(*failHisto) );
//	RooPlot* failFrame = fx.frame(Title("fail TH1 with Poisson error bars"));
//	failData.plotOn(failFrame);
//	RooGaussian fgauss("fgauss", "fgaussian PDF", fx, fmean, fsigma) ;
//
//	fgauss.fitTo (failData );
//	fgauss.plotOn(failFrame);


	// drawing things
 	TLatex *lat = new TLatex();
 	lat->SetNDC(kTRUE);
 	lat->SetTextColor(kBlack);
 	lat->SetTextSize(0.05);

 	TString indText;
 	if     (flag==0) indText = "ID over all";
 	else if(flag==1) indText = "IP over ID";
 	else if(flag==2) indText = "ISO over ID+IP";
 
 	TCanvas* c = new TCanvas("foo", "bar", 800, 675);
 	c->Divide(1,2);
 	c->cd(1);
	passFrame->Draw();
 	lat->DrawLatex(0.10,0.92, indText+" passing");
 	lat->DrawLatex(0.40,0.92, fIsMu ? getPtEtaFromIndexMu(bin) : getPtEtaFromIndexEl(bin));

// 	c->cd(2);
//	failFrame->Draw();
// 	lat->DrawLatex(0.10,0.92, indText+" failing");
// 	lat->DrawLatex(0.40,0.92, fIsMu ? getPtEtaFromIndexMu(bin) : getPtEtaFromIndexEl(bin));


 	// save as root and pdf
 	c->Print(Form("passAndFail_bin%d.pdf" ,bin));
 	c->Print(Form("passAndFail_bin%d.root",bin));

	delete passFrame, c, lat;
	//delete failFrame;

}

// void TnP::simFitPassFail(TH1F* passHisto , TH1F* failHisto, int flag, int bin){
// 
// 	// linear + gauss
// 	TF1 * fP = new TF1("fP","[0] + [1]*x + [2] * 1/(sqrt(2*TMath::Pi())*[4]) * exp(-0.5*((x-[3])/[4])**2)", 60, 120);
// 	// exponential + gauss
// 	// TF1 * fP = new TF1("fP","[0]*exp([5] + [1]*x) + [2] * 1/(sqrt(2*TMath::Pi())*[4]) * exp(-0.5*((x-[3])/[4])**2)", 60, 120);
// 	fP->SetNpx(1000);
// 	fP->SetParameter(3,90);
// 	fP->SetParameter(4, 5);
// 
// 	fP->SetLineColor(4);
// 	passHisto->Fit("fP", "eq");
// 	float npass    = fP->GetParameter(2);
// 	float npassErr = fP->GetParError(2);
// 
// 	// fix the mean and the sigma of the gaussian
// 	fP->FixParameter(3, fP->GetParameter(3));
// 	fP->FixParameter(4, fP->GetParameter(4));
// 	
// 
// 	fP->SetLineColor(2);
// 	failHisto->Fit("fP", "q");
// 	float nfail    = fP->GetParameter(2);
// 	float nfailErr = fP->GetParError(2);
// 
// 	float eff = (npass+nfail) > 0 ? npass/(npass+nfail) : 0.;
// 	float effErr = (npass+nfail) > 0 ? fabs(npass/(npass+nfail)*(npassErr/npass - nfailErr/nfail)) : 0.;
// 
// 	TLatex *lat = new TLatex();
// 	lat->SetNDC(kTRUE);
// 	lat->SetTextColor(kBlack);
// 	lat->SetTextSize(0.05);
// 
// 	TString indText;
// 	if     (flag==0) indText = "ID over all";
// 	else if(flag==1) indText = "IP over ID";
// 	else if(flag==2) indText = "ISO over ID+IP";
// 
// 	TCanvas* c = new TCanvas("foo", "bar", 800, 675);
// 	c->Divide(1,2);
// 	c->cd(1);
// 	passHisto->Draw("pe");
// 	setHistoStyle(passHisto);
// 	passHisto->GetXaxis()->SetRangeUser(60., 120.);
// 	lat->DrawLatex(0.10,0.92, indText+" passing");
// 	lat->DrawLatex(0.40,0.92, fIsMu ? getPtEtaFromIndexMu(bin) : getPtEtaFromIndexEl(bin));
// 
// 	c->cd(2);
// 	failHisto->Draw("pe");
// 	setHistoStyle(failHisto);
// 	failHisto->GetXaxis()->SetRangeUser(60., 120.);
// 	lat->DrawLatex(0.10,0.92, indText+" failing");
// 	lat->DrawLatex(0.40,0.92, fIsMu ? getPtEtaFromIndexMu(bin) : getPtEtaFromIndexEl(bin));
// 
// 	// save as root and pdf
// 	c->Print(Form("passAndFail_bin%d.pdf" ,bin));
// 	c->Print(Form("passAndFail_bin%d.root",bin));
// 
// 	delete c;
// 
// 	// store the efficiencies:
// 	effs[bin].bin = bin;
// 	if(flag == 0) {
// 		effs[bin].idEff = eff;
// 		effs[bin].idEffErr = effErr;
// 	}
// 	if(flag == 1) {
// 		effs[bin].ipEff = eff;
// 		effs[bin].ipEffErr = effErr;
// 	}
// 	if(flag == 2) {
// 		effs[bin].isoEff = eff;
// 		effs[bin].isoEffErr = effErr;
// 	}
// 
// }

void TnP::fillHistos(){

	// load the tree

	TFile *pFile = TFile::Open(fInputFile);
	TTree *tree; getObjectSafe(pFile, "probeTree", tree);
	float tagPt, tagEta, tagPhi, tagIsoRel, tagD0, tagDz;
	float tagPassID;
	float probePt, probeEta, probePhi, probeIsoRel, probeD0, probeDz;
	float probePassID;
	float mll, deltaR, nvtx, pfmet;// , rho, ht, nTrueInt;
	// tag variables
	tree->SetBranchAddress("tagPt"       , &tagPt       );
	tree->SetBranchAddress("tagEta"      , &tagEta      );
	tree->SetBranchAddress("tagPhi"      , &tagPhi      );
	tree->SetBranchAddress("tagIsoRel"   , &tagIsoRel   );
	tree->SetBranchAddress("tagD0"       , &tagD0       );
	tree->SetBranchAddress("tagDz"       , &tagDz       );
	tree->SetBranchAddress("tagPassID"   , &tagPassID   );
	// probe variables
	tree->SetBranchAddress("probePt"     , &probePt     );
	tree->SetBranchAddress("probeEta"    , &probeEta    );
	tree->SetBranchAddress("probePhi"    , &probePhi    );
	tree->SetBranchAddress("probeIsoRel" , &probeIsoRel );
	tree->SetBranchAddress("probeD0"     , &probeD0     );
	tree->SetBranchAddress("probeDz"     , &probeDz     );
	tree->SetBranchAddress("probePassID" , &probePassID );
	// other variables
	tree->SetBranchAddress("mass2L"      , &mll         );
	tree->SetBranchAddress("deltaR"      , &deltaR      );
	tree->SetBranchAddress("nvtx"        , &nvtx        );
	tree->SetBranchAddress("pfmet"       , &pfmet       );

	long nentries = tree->GetEntries();
	cout << "looping over " << nentries << " events" << endl;

	int index = -1;	

	for( int i = 0; i < nentries; i++ ){
		tree->GetEntry(i);
		if(!passesAll(tagIsoRel, tagD0, tagDz, tagPassID)) continue; // make sure the tag passes everything
		index = fIsMu ? getPtEtaIndexMu(probePt, probeEta) : getPtEtaIndexEl(probePt, probeEta);
		if(!fIsMu and !checkElEta(probeEta)) continue;
		if(probePassID == 0){                               // probe fails ID
			fFailIDoverAll[index]   ->Fill(mll);            // fill the histo
			fFailIDoverAll[fnBins-1]->Fill(mll);            // fill the inclusive histo
		}
		else {                                              // probe passes ID
			fPassIDoverAll[index]   ->Fill(mll);            // fill the histo
			fPassIDoverAll[fnBins-1]->Fill(mll);            // fill the inclusive histo
			if(!passesIP(probeD0, probeDz)){                // probe fails IP but passes ID
				fFailIPoverID[index]   ->Fill(mll);         // fill the histo
				fFailIPoverID[fnBins-1]->Fill(mll);         // fill the inclusive histo
			}
			else {                                          // probe passes ID+IP
				fPassIPoverID[index]   ->Fill(mll);         // fill the histo
				fPassIPoverID[fnBins-1]->Fill(mll);         // fill the inclusive histo
				if(!passesIso(probeIsoRel)){                // probe fails ISO but passes ID+IP
					fFailISOoverIDIP[index]   ->Fill(mll);  // fill the histo
					fFailISOoverIDIP[fnBins-1]->Fill(mll);  // fill the inclusive histo
				}
				else{                                       // probe passes ID+IP+ISO
					fPassISOoverIDIP[index]   ->Fill(mll);  // fill the histo
					fPassISOoverIDIP[fnBins-1]->Fill(mll);  // fill the inclusive histo
				}
			}
		}
	}
	pFile->Close();

}
bool TnP::passesAll(float iso, float d0, float dz, int passid){
	if(iso       > fIsoCut) return false;
	if(fabs(d0)  > fD0Cut ) return false;
	if(fabs(dz)  > 0.10  ) return false;
	if(passid   == 0     ) return false;
	return true;
}
bool TnP::passesID(int passid){
	if(passid   == 0) return false;
	return true;
}
bool TnP::passesIP(float d0, float dz){
	if(fabs(d0)  > fD0Cut) return false;
	if(fabs(dz)  > 0.10 ) return false;
	return true;
}
bool TnP::passesIso(float iso){
	if(iso       > fIsoCut) return false;
	return true;
}
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________

void TnP::writeHistos(){
	
	TFile *pFile = TFile::Open(fInputFile, "UPDATE");
	pFile->cd();
	for(int i=0; i<fnBins; ++i){
		fPassIDoverAll   [i]->Write(fPassIDoverAll   [i]->GetName(), TObject::kOverwrite);
		fFailIDoverAll   [i]->Write(fFailIDoverAll   [i]->GetName(), TObject::kOverwrite);

		fPassIPoverID    [i]->Write(fPassIPoverID    [i]->GetName(), TObject::kOverwrite);
		fFailIPoverID    [i]->Write(fFailIPoverID    [i]->GetName(), TObject::kOverwrite);

		fPassISOoverIDIP [i]->Write(fPassISOoverIDIP [i]->GetName(), TObject::kOverwrite);
		fFailISOoverIDIP [i]->Write(fFailISOoverIDIP [i]->GetName(), TObject::kOverwrite);
	}
	pFile->Close();
}
void TnP::bookHistos(){
}
void TnP::deleteHistos(){
}
int  TnP::readHistos(TString filename){
	std::cout << "this is the filename of which to read the histos: " << filename << std::endl;
	return 0;
}

// =====================================================
// MUON FUNCTIONS
// =====================================================

int TnP::getPtEtaIndexMu(float pt, float eta){
	float aeta = fabs(eta);
	if(pt > 10. and pt <= 15.){
		if(aeta <= 1.2) return 0;
		if(aeta >  1.2) return 6;
	}
	else if(pt > 15. and pt <= 20.){
		if(aeta <= 1.2) return 1;
		if(aeta >  1.2) return 7;
	}
	else if(pt > 20. and pt <= 30.){
		if(aeta <= 1.2) return 2;
		if(aeta >  1.2) return 8;
	}
	else if(pt > 30. and pt <= 40.){
		if(aeta <= 1.2) return 3;
		if(aeta >  1.2) return 9;
	}
	else if(pt > 40. and pt <= 50.){
		if(aeta <= 1.2) return 4;
		if(aeta >  1.2) return 10;
	}
	else if(pt > 50.){
		if(aeta <= 1.2) return 5;
		if(aeta >  1.2) return 11;
	}
	else{
		cout << " ======== ERROR ================= " << endl;
		cout << " something went wrong with the index return for your muon... " << endl;
		cout << " exiting... " << endl;
		exit(-1);
	}
	return -1;
}
TString TnP::getPtEtaFromIndexMu(int i){
	TString str;
	switch (i) {
		case  0: str = "10 < pT < 15 GeV  ; |eta| <= 1.2"; break;
		case  6: str = "10 < pT < 15 GeV  ; |eta| >  1.2"; break;
		case  1: str = "15 < pT < 20 GeV  ; |eta| <= 1.2"; break;
		case  7: str = "15 < pT < 20 GeV  ; |eta| >  1.2"; break;
		case  2: str = "20 < pT < 30 GeV  ; |eta| <= 1.2"; break;
		case  8: str = "20 < pT < 30 GeV  ; |eta| >  1.2"; break;
		case  3: str = "30 < pT < 40 GeV  ; |eta| <= 1.2"; break;
		case  9: str = "30 < pT < 40 GeV  ; |eta| >  1.2"; break;
		case  4: str = "40 < pT < 50 GeV  ; |eta| <= 1.2"; break;
		case 10: str = "40 < pT < 50 GeV  ; |eta| >  1.2"; break;
		case  5: str = "     pT > 50 GeV  ; |eta| <= 1.2"; break;
		case 11: str = "     pT > 50 GeV  ; |eta| >  1.2"; break;
		case 12: str = "    inclusive                   "; break;
	}
	return str;
}
void TnP::printMuTable(){
	cout << Form("========================================================") << endl;
	cout << Form("========================================================") << endl;
	for(int i=0; i<fnBins; ++i){
		cout << Form("At bin %2d: "+getPtEtaFromIndexMu(i), i+1) << endl;
		cout << Form("========================") << endl;
		cout << Form("ID-eff : %.3f +- %.3f", effs[i].idEff , effs[i].idEffErr ) << endl;
		cout << Form("IP-eff : %.3f +- %.3f", effs[i].ipEff , effs[i].ipEffErr ) << endl;
		cout << Form("ISO-eff: %.3f +- %.3f", effs[i].isoEff, effs[i].isoEffErr) << endl;
		cout << Form("=============================================") << endl;
	}
}

// =====================================================
// ELECTRON FUNCTIONS
// =====================================================
bool TnP::checkElEta(float eta){
	float aeta = fabs(eta);
	if(aeta > 2.4) return false;
	if(aeta > 1.4442 and aeta < 1.566) return false;
	return true;
}

int TnP::getPtEtaIndexEl(float pt, float eta){
	float aeta = fabs(eta);
	// cout << Form("pt: %.2f   eta: %.2f", pt, aeta) << endl;
	if     (pt >= 10. and pt < 15.){
		if( 0.0   <= aeta and aeta <  0.8   ) return  0;
		if( 0.8   <= aeta and aeta <  1.4442) return  6;
		if( 1.566 <= aeta and aeta <  2.0   ) return 12;
		if( 2.0   <= aeta and aeta <  2.4   ) return 18;
	}
	else if(pt >= 15. and pt < 20.){
		if( 0.0   <= aeta and aeta <  0.8   ) return  1;
		if( 0.8   <= aeta and aeta <  1.4442) return  7;
		if( 1.566 <= aeta and aeta <  2.0   ) return 13;
		if( 2.0   <= aeta and aeta <  2.4   ) return 19;
	}
	else if(pt >= 20. and pt < 30.){
		if( 0.0   <= aeta and aeta <  0.8   ) return  2;
		if( 0.8   <= aeta and aeta <  1.4442) return  8;
		if( 1.566 <= aeta and aeta <  2.0   ) return 14;
		if( 2.0   <= aeta and aeta <  2.4   ) return 20;
	}
	else if(pt >= 30. and pt < 40.){
		if( 0.0   <= aeta and aeta <  0.8   ) return  3;
		if( 0.8   <= aeta and aeta <  1.4442) return  9;
		if( 1.566 <= aeta and aeta <  2.0   ) return 15;
		if( 2.0   <= aeta and aeta <  2.4   ) return 21;
	}
	else if(pt >= 40. and pt < 50.){
		if( 0.0   <= aeta and aeta <  0.8   ) return  4;
		if( 0.8   <= aeta and aeta <  1.4442) return 10;
		if( 1.566 <= aeta and aeta <  2.0   ) return 16;
		if( 2.0   <= aeta and aeta <  2.4   ) return 22;
	}
	else if(pt >= 50.){
		if( 0.0   <= aeta and aeta <  0.8   ) return  5;
		if( 0.8   <= aeta and aeta <  1.4442) return 11;
		if( 1.566 <= aeta and aeta <  2.0   ) return 17;
		if( 2.0   <= aeta and aeta <  2.4   ) return 23;
	}
	else{
		cout << " ======== ERROR ================= " << endl;
		cout << " something went wrong with the index return for your muon... " << endl;
		cout << " exiting... " << endl;
		exit(-1);
	}
	return -1;
}
TString TnP::getPtEtaFromIndexEl(int i){
	TString str;
	switch (i) {
		case  0: str = "10 < pT < 15 GeV  ; |eta| in [0.0    , 0.8   ]"; break;
		case  6: str = "10 < pT < 15 GeV  ; |eta| in [0.8    , 1.4442]"; break;
		case 12: str = "10 < pT < 15 GeV  ; |eta| in [1.556  , 2.0   ]"; break;
		case 18: str = "10 < pT < 15 GeV  ; |eta| in [2.0    , 2.4   ]"; break;

		case  1: str = "15 < pT < 20 GeV  ; |eta| in [0.0    , 0.8   ]"; break;
		case  7: str = "15 < pT < 20 GeV  ; |eta| in [0.8    , 1.4442]"; break;
		case 13: str = "15 < pT < 20 GeV  ; |eta| in [1.556  , 2.0   ]"; break;
		case 19: str = "15 < pT < 20 GeV  ; |eta| in [2.0    , 2.4   ]"; break;

		case  2: str = "20 < pT < 30 GeV  ; |eta| in [0.0    , 0.8   ]"; break;
		case  8: str = "20 < pT < 30 GeV  ; |eta| in [0.8    , 1.4442]"; break;
		case 14: str = "20 < pT < 30 GeV  ; |eta| in [1.556  , 2.0   ]"; break;
		case 20: str = "20 < pT < 30 GeV  ; |eta| in [2.0    , 2.4   ]"; break;

		case  3: str = "30 < pT < 40 GeV  ; |eta| in [0.0    , 0.8   ]"; break;
		case  9: str = "30 < pT < 40 GeV  ; |eta| in [0.8    , 1.4442]"; break;
		case 15: str = "30 < pT < 40 GeV  ; |eta| in [1.556  , 2.0   ]"; break;
		case 21: str = "30 < pT < 40 GeV  ; |eta| in [2.0    , 2.4   ]"; break;

		case  4: str = "40 < pT < 50 GeV  ; |eta| in [0.0    , 0.8   ]"; break;
		case 10: str = "40 < pT < 50 GeV  ; |eta| in [0.8    , 1.4442]"; break;
		case 16: str = "40 < pT < 50 GeV  ; |eta| in [1.556  , 2.0   ]"; break;
		case 22: str = "40 < pT < 50 GeV  ; |eta| in [2.0    , 2.4   ]"; break;

		case  5: str = "     pT > 50 GeV  ; |eta| in [0.0    , 0.8   ]"; break;
		case 11: str = "     pT > 50 GeV  ; |eta| in [0.8    , 1.4442]"; break;
		case 17: str = "     pT > 50 GeV  ; |eta| in [1.556  , 2.0   ]"; break;
		case 23: str = "     pT > 50 GeV  ; |eta| in [2.0    , 2.4   ]"; break;

		case 24: str = "               inclusive                      "; break;
	}
	return str;
}
void TnP::printElTable(){
	cout << Form("========================================================") << endl;
	cout << Form("========================================================") << endl;
	for(int i=0; i<fnBins; ++i){
		cout << Form("At bin %2d: "+getPtEtaFromIndexEl(i), i+1) << endl;
		cout << Form("========================") << endl;
		cout << Form("ID-eff : %.3f +- %.3f", effs[i].idEff , effs[i].idEffErr ) << endl;
		cout << Form("IP-eff : %.3f +- %.3f", effs[i].ipEff , effs[i].ipEffErr ) << endl;
		cout << Form("ISO-eff: %.3f +- %.3f", effs[i].isoEff, effs[i].isoEffErr) << endl;
		cout << Form("=============================================") << endl;
	}
}

// ------------------------------------------
// --  STYLE STUFF  -------------------------
// ------------------------------------------

void TnP::setHistoStyle(TH1F* histo){
	histo->SetMarkerStyle(20);
	histo->SetMarkerSize(0.9);
	histo->SetLineWidth(0.8);
	gStyle->SetOptStat(1111111);
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);
}

