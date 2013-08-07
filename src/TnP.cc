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


// stuff for fitting
#include "TF1.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <iostream>
#include <iomanip>
#include <time.h> // access to date/time

#include <typeinfo>


using namespace std;

//____________________________________________________________________________
TnP::TnP(TString inputfile, bool createHistos){
	fCreateHistos = createHistos;
	fInputFile = inputfile;

	fOutputFileName = "tnp_output.root";

	fMuIso = 0.09;
	fMuD0  = 0.005;


	if(fCreateHistos){
		for(int i=0; i<fnBinsMu; ++i){
			fMuPassIDoverAll   [i] = new TH1D(Form("MuPassIDoverAll_%d"  ,i), Form("MuPassIDoverAll_%d"  ,i), 120, 0., 120.);
			fMuPassIPoverID    [i] = new TH1D(Form("MuPassIPoverID_%d"   ,i), Form("MuPassIPoverID_%d"   ,i), 120, 0., 120.);
			fMuPassISOoverIDIP [i] = new TH1D(Form("MuPassISOoverIDIP_%d",i), Form("MuPassISOoverIDIP_%d",i), 120, 0., 120.);
                                                                                                  
			fMuFailIDoverAll   [i] = new TH1D(Form("MuFailIDoverAll_%d"  ,i), Form("MuFailIDoverAll_%d"  ,i), 120, 0., 120.);
			fMuFailIPoverID    [i] = new TH1D(Form("MuFailIPoverID_%d"   ,i), Form("MuFailIPoverID_%d"   ,i), 120, 0., 120.);
			fMuFailISOoverIDIP [i] = new TH1D(Form("MuFailISOoverIDIP_%d",i), Form("MuFailISOoverIDIP_%d",i), 120, 0., 120.);
		}
	}
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

		for(int i=0; i<fnBinsMu; ++i){
			if(fVerbose > 2) inFile->ls();
			getObjectSafe(inFile, Form("MuPassIDoverAll_%d"  ,i), fMuPassIDoverAll   [i]);
			getObjectSafe(inFile, Form("MuPassIPoverID_%d"   ,i), fMuPassIPoverID    [i]);
			getObjectSafe(inFile, Form("MuPassISOoverIDIP_%d",i), fMuPassISOoverIDIP [i]);

			getObjectSafe(inFile, Form("MuFailIDoverAll_%d"  ,i), fMuFailIDoverAll   [i]);
			getObjectSafe(inFile, Form("MuFailIPoverID_%d"   ,i), fMuFailIPoverID    [i]);
			getObjectSafe(inFile, Form("MuFailISOoverIDIP_%d",i), fMuFailISOoverIDIP [i]);
		}
		
	}
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

	for(int i=0; i<fnBinsMu; ++i){
		cout << "======================================" << endl;
		cout << "at bin: " << i << endl;
		cout << "======================================" << endl;
		cout << "ID:" << endl << "===" << endl;
		simFitPassFail(fMuPassIDoverAll[i], fMuFailIDoverAll[i], 0, i);
		cout << "IP:" << endl << "===" << endl;
		simFitPassFail(fMuPassIPoverID[i], fMuFailIPoverID[i], 0, i);
		cout << "ISO:" << endl << "====" << endl;
		simFitPassFail(fMuPassISOoverIDIP[i], fMuFailISOoverIDIP[i], 0, i);
	}

}

void TnP::simFitPassFail(TH1D* passHisto , TH1D* failHisto, int flag, int bin){

	TF1 * fP = new TF1("fP","[0] + [1]*x + [2] * 1/(sqrt(2*TMath::Pi())*[4]) * exp(-0.5*((x-[3])/[4])**2)", 0, 120);
	fP->SetParameter(3,90);
	fP->SetParameter(4, 5);

	passHisto->Fit("fP", "eq");
	float npass    = fP->GetParameter(2);
	float npassErr = fP->GetParError(2);

	// fix the mean and the sigma of the gaussian
	fP->FixParameter(3, fP->GetParameter(3));
	fP->FixParameter(4, fP->GetParameter(4));
	

	failHisto->Fit("fP", "q");
	float nfail    = fP->GetParameter(2);
	float nfailErr = fP->GetParError(2);

	float eff = (npass+nfail) > 0 ? npass/(npass+nfail) : 0.;
	float effErr = (npass+nfail) > 0 ? fabs(npass/(npass+nfail)*(npassErr/npass - nfailErr/nfail)) : 0.;

	// std::cout << Form("npass: %.2f  +- %.2f nfail: %.2f +- %.2f", npass, npassErr, nfail, nfailErr) << std::endl;
	// std::cout << Form("efficiency: %.3f +- %.3f", eff, effErr)   << std::endl;
	
	muEffs[bin].bin = bin;
	if(flag == 0) {
		muEffs[bin].idEff = eff;
		muEffs[bin].idEffErr = effErr;
	}
	if(flag == 1) {
		muEffs[bin].ipEff = eff;
		muEffs[bin].ipEffErr = effErr;
	}
	if(flag == 2) {
		muEffs[bin].isoEff = eff;
		muEffs[bin].isoEffErr = effErr;
	}

}

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
		if(!passesAllMu(tagIsoRel, tagD0, tagDz, tagPassID)) continue; // make sure the tag passes everything
		index = getPtEtaIndexMu(probePt, probeEta);
		if(probePassID == 0){                                   // probe fails ID
			fMuFailIDoverAll[index]     ->Fill(mll);            // fill the histo
			fMuFailIDoverAll[fnBinsMu-1]->Fill(mll);            // fill the inclusive histo
		}
		else {                                                  // probe passes ID
			fMuPassIDoverAll[index]     ->Fill(mll);            // fill the histo
			fMuPassIDoverAll[fnBinsMu-1]->Fill(mll);            // fill the inclusive histo
			if(!passesIPMu(probeD0, probeDz)){                  // probe fails IP but passes ID
				fMuFailIPoverID[index]     ->Fill(mll);         // fill the histo
				fMuFailIPoverID[fnBinsMu-1]->Fill(mll);         // fill the inclusive histo
			}
			else {                                              // probe passes ID+IP
				fMuPassIPoverID[index]     ->Fill(mll);         // fill the histo
				fMuPassIPoverID[fnBinsMu-1]->Fill(mll);         // fill the inclusive histo
				if(!passesIsoMu(probeIsoRel)){                  // probe fails ISO but passes ID+IP
					fMuFailISOoverIDIP[index]     ->Fill(mll);  // fill the histo
					fMuFailISOoverIDIP[fnBinsMu-1]->Fill(mll);  // fill the inclusive histo
				}
				else{                                           // probe passes ID+IP+ISO
					fMuPassISOoverIDIP[index]     ->Fill(mll);  // fill the histo
					fMuPassISOoverIDIP[fnBinsMu-1]->Fill(mll);  // fill the inclusive histo
				}
			}
		}
	}
	pFile->Close();

}
void TnP::writeHistos(){
	
	TFile *pFile = TFile::Open(fInputFile, "UPDATE");
	pFile->cd();
	for(int i=0; i<fnBinsMu; ++i){
		fMuPassIDoverAll   [i]->Write(fMuPassIDoverAll   [i]->GetName(), TObject::kOverwrite);
		fMuFailIDoverAll   [i]->Write(fMuFailIDoverAll   [i]->GetName(), TObject::kOverwrite);

		fMuPassIPoverID    [i]->Write(fMuPassIPoverID    [i]->GetName(), TObject::kOverwrite);
		fMuFailIPoverID    [i]->Write(fMuFailIPoverID    [i]->GetName(), TObject::kOverwrite);

		fMuPassISOoverIDIP [i]->Write(fMuPassISOoverIDIP [i]->GetName(), TObject::kOverwrite);
		fMuFailISOoverIDIP [i]->Write(fMuFailISOoverIDIP [i]->GetName(), TObject::kOverwrite);
	}
	pFile->Close();
}

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
// int TnP::getPtEtaIndexEl(float pt, float eta){
// 	return 1;
// }

void TnP::bookHistos(){
}
void TnP::deleteHistos(){
}
int  TnP::readHistos(TString filename){
	std::cout << "this is the filename of which to read the histos: " << filename << std::endl;
	return 0;
}

//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________

bool TnP::passesAllMu(float iso, float d0, float dz, int passid){
	if(iso       > fMuIso) return false;
	if(fabs(d0)  > fMuD0 ) return false;
	if(fabs(dz)  > 0.10  ) return false;
	if(passid   == 0     ) return false;
	return true;
}
bool TnP::passesIDMu(int passid){
	if(passid   == 0) return false;
	return true;
}
bool TnP::passesIPMu(float d0, float dz){
	if(fabs(d0)  > fMuD0) return false;
	if(fabs(dz)  > 0.10 ) return false;
	return true;
}
bool TnP::passesIsoMu(float iso){
	if(iso       > fMuIso) return false;
	return true;
}
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________

