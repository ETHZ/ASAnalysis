#include "helper/Utilities.hh"
#include "TnPAnalysis.hh"

using namespace std;

TnPAnalysis::TnPAnalysis(TreeReader *tr) : UserAnalysisBase(tr){
    Util::SetStyle();
}

TnPAnalysis::~TnPAnalysis(){
}

void TnPAnalysis::Begin(){
	// Define the output file of histograms
	//const char* filename = "histos.root";
	//fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	fHistFile     = new TFile(fOutputFileName, "RECREATE");
	
	// Define the histograms
	fHCounter     = new TH1D("counter"   , "counter"   ,   1,     0,     1);

	fHPtHlt1      = new TH1D("ptHlt1"    , "ptHlt1"    , 100,    0.,  100.);
	fHEtaHlt1     = new TH1D("etaHlt1"   , "etaHlt1"   ,  50,    -4,    4.);
	fHPhiHlt1     = new TH1D("phiHlt1"   , "phiHlt1"   ,  50, -3.14,  3.14);

	fHPtHlt2      = new TH1D("ptHlt2"    , "ptHlt2"    , 100,    0.,  100.);
	fHEtaHlt2     = new TH1D("etaHlt2"   , "etaHlt2"   ,  50,    -4,    4.);
	fHPhiHlt2     = new TH1D("phiHlt2"   , "phiHlt2"   ,  50, -3.14,  3.14);

	fHDeltaR      = new TH1D("deltaR"    , "deltaR"    , 100,    0., 0.005);
	fHPtTag       = new TH1D("ptTag"     , "ptTag"     , 100,    0.,  100.);
	fHMass2LAll   = new TH1D("mass2LAll" , "mass2LAll" , 120,    0.,  120.);
	fHMass2LPass  = new TH1D("mass2LPass", "mass2LPass", 120,    0.,  120.);
	fHMass2LFail  = new TH1D("mass2LFail", "mass2LFail", 120,    0.,  120.);

	BookTree();
}

void TnPAnalysis::SetIsMu(bool isMu){
	fIsMu = isMu;
}

void TnPAnalysis::Analyze(){
	fHCounter->Fill(0.5);
	// Some event selection
	
	if(fIsMu ? !fTR->HLTObjectPt0.size() : !fTR->HLTObjectPt1.size()) {
		// cout << "warning: hltObjPt is empty" << endl;
		return;
	}
	
	LookForTagAndProbe();
	
}

void TnPAnalysis::LookForTagAndProbe(){

	float massLep = fIsMu ? 0.10566 : 0.000511;

	// ----- look for tag
	TLorentzVector lepHLT;
	lepHLT.SetPtEtaPhiM(fIsMu ? fTR->HLTObjectPt0 [0] : fTR->HLTObjectPt1 [0], 
	  	                fIsMu ? fTR->HLTObjectEta0[0] : fTR->HLTObjectEta1[0], 
	  	                fIsMu ? fTR->HLTObjectPhi0[0] : fTR->HLTObjectPhi1[0], 
	  	                massLep);
	
	
	TLorentzVector lepTag;
	float deltaR(9999.); 
	int   indexLepTag(-1);
	bool  foundTag(false);
	int   loopMax = fIsMu ? fTR->NMus : fTR->NEles;
	for(int i=0; i< loopMax; i++){
		//cout << "mu pt: " << fTR->MuPt[i] << endl;
		lepTag.SetPtEtaPhiM(fIsMu ? fTR->MuPt [i] : fTR->ElPt [i],
				            fIsMu ? fTR->MuEta[i] : fTR->ElEta[i],
				            fIsMu ? fTR->MuPhi[i] : fTR->ElPhi[i],
				            massLep);
		float tmpDeltaR = lepHLT.DeltaR(lepTag);
		float deltaP    = fabs(lepHLT.Pt() - lepTag.Pt())/lepHLT.Pt();
		float lepIso    = fIsMu ? MuPFIso(i)    : ElPFIso(i);
		bool  passID    = fIsMu ? MuPassesPOGTightID(i) : ElPassesPOGMediumWP(i);
		float isoThresh = fIsMu ? 0.1 : 0.09;
		if( tmpDeltaR < deltaR && deltaP < 0.05 && passID  && lepIso < isoThresh ){
			deltaR = tmpDeltaR;
			indexLepTag = i;
			foundTag = true;
			break;
		}
	}

	if(foundTag){
		fHDeltaR->Fill(deltaR);
		fHPtTag ->Fill(fIsMu ? fTR->MuPt[indexLepTag] : fTR->ElPt[indexLepTag]);
	}

	// the probe
	TLorentzVector lepProbe;
	bool foundProbe(false); int indexLepProbe(-1);
	if(foundTag){
		for(int i=0; i< loopMax; i++){
			if(i==indexLepTag) continue;
			lepProbe.SetPtEtaPhiM(fIsMu ? fTR->MuPt [i] : fTR->ElPt [i],
			  	                  fIsMu ? fTR->MuEta[i] : fTR->ElEta[i],
			  	                  fIsMu ? fTR->MuPhi[i] : fTR->ElPhi[i],
			  	                  massLep);
			float tmpDeltaR = lepTag.DeltaR(lepProbe);
			if(tmpDeltaR     < 0.1) continue;
			if(lepProbe.Pt() < 10.) continue;
			foundProbe    = true;
			indexLepProbe = i;
			break;
		}
	}
	// ---
	
	// filling tree after tag & probe are found
	if(foundTag && foundProbe){
		float mass2L = (lepTag+lepProbe).M();
		fHMass2LAll->Fill(mass2L);
		if(fIsMu ? MuPassesPOGTightID(indexLepProbe) : ElPassesPOGMediumWP(indexLepProbe))
			fHMass2LPass->Fill(mass2L);
		else
			fHMass2LFail->Fill(mass2L);
		FillAnalysisTree(indexLepTag, indexLepProbe, massLep);
	}
	
	
	
	// if(fTR->HLTObjectPt0.size()>=1){
	// 	fIsMu ? fHPtHlt1 ->Fill(fTR->HLTObjectPt0 [0]) : fHPtHlt1 ->Fill(fTR->HLTObjectPt1 [0]);
	// 	fIsMu ? fHEtaHlt1->Fill(fTR->HLTObjectEta0[0]) : fHEtaHlt1->Fill(fTR->HLTObjectEta1[0]);
	// 	fIsMu ? fHPhiHlt1->Fill(fTR->HLTObjectPhi0[0]) : fHPhiHlt1->Fill(fTR->HLTObjectPhi1[0]);
	// }
	// if(fTR->HLTObjectPt1.size()>=1){
	// 	fIsMu ? fHPtHlt2 ->Fill(fTR->HLTObjectPt1 [1]);
	// 	fIsMu ? fHEtaHlt2->Fill(fTR->HLTObjectEta1[1]);
	// 	fIsMu ? fHPhiHlt2->Fill(fTR->HLTObjectPhi1[1]);
	// }
	
}

void TnPAnalysis::End(){
	fHistFile->cd();
	fHCounter->Write();	
	
	fHPtHlt1 ->Write();
	fHEtaHlt1 ->Write();
	fHPhiHlt1 ->Write();
	
	fHPtHlt2 ->Write();
	fHEtaHlt2 ->Write();
	fHPhiHlt2 ->Write();
	
	fHDeltaR->Write();
	fHPtTag->Write();
	fHMass2LAll->Write();
	fHMass2LPass->Write();
	fHMass2LFail->Write();
	
	fAnalysisTree->Write();
	
	fHistFile->Close();
	//fOutputFile->cd();
	//fOutputFile->Close();

}




void TnPAnalysis::BookTree(){
	//cout << "here 1" << endl;
	//cout << "fOutputfile: " << fOutputFile << endl;
	//fOutputFile->cd();
	//cout << "here 2" << endl;
	fAnalysisTree = new TTree("probeTree", "probeTree");
	
	// run/sample properties
	fAnalysisTree->Branch("Run",              &fTRunNumber,         "Run/I");
	fAnalysisTree->Branch("Event",            &fTEventNumber,       "Event/I");
	fAnalysisTree->Branch("LumiSec",          &fTLumiSection,       "LumiSec/I");
	
	fAnalysisTree->Branch("isData",         &fTisData,     "isData/I");
	fAnalysisTree->Branch("isMuEvent",         &fTisMuEvent,     "isMuEvent/I");
	
	fAnalysisTree->Branch("tagPt",            &fTagPt,         "tagPt/F");
	fAnalysisTree->Branch("tagEta",           &fTagEta,        "tagEta/F");
	fAnalysisTree->Branch("tagPhi",           &fTagPhi,        "tagPhi/F");
	fAnalysisTree->Branch("tagIsoRel",        &fTagIsoRel,     "tagIsoRel/F");
	fAnalysisTree->Branch("tagD0",            &fTagD0,         "tagD0/F");
	fAnalysisTree->Branch("tagDz",            &fTagDz,         "tagDz/F");
	fAnalysisTree->Branch("tagPassID",        &fTagPassID,     "tagPassID/F");
	//
	fAnalysisTree->Branch("probePt",            &fProbePt,         "probePt/F");
	fAnalysisTree->Branch("probeEta",           &fProbeEta,        "probeEta/F");
	fAnalysisTree->Branch("probePhi",           &fProbePhi,        "probePhi/F");
	fAnalysisTree->Branch("probeIsoRel",        &fProbeIsoRel,     "probeIsoRel/F");
	fAnalysisTree->Branch("probeD0",            &fProbeD0,         "probeD0/F");
	fAnalysisTree->Branch("probeDz",            &fProbeDz,         "probeDz/F");
	fAnalysisTree->Branch("probePassID",        &fProbePassID,     "probePassID/F");
	
	fAnalysisTree->Branch("mass2L",        &fMass2L,     "mass2L/F");
	fAnalysisTree->Branch("deltaR",        &fDeltaR,     "deltaR/F");
	fAnalysisTree->Branch("nvtx",          &fNvtx,       "nvtx/F");
	fAnalysisTree->Branch("rho",           &fRho,        "rho/F");
	fAnalysisTree->Branch("pfmet",         &fPfMet,      "pfmet/F");
	fAnalysisTree->Branch("ht",            &fHt,         "ht/F");
	fAnalysisTree->Branch("nTrueInt",      &fNTrueInt,   "nTrueInt/I");

}

void TnPAnalysis::FillAnalysisTree(int indexLepTag, int indexLepProbe, float massLep){

	// Event and run info
	fTRunNumber   = fTR->Run;
	fTEventNumber = fTR->Event;
	fTLumiSection = fTR->LumiSection;
	
	fTisData = fIsData ? 1 : 0;
	
	
	TLorentzVector lepTag;
	TLorentzVector lepProbe;
	
	
	lepTag  .SetPtEtaPhiM(fIsMu ? fTR->MuPt [indexLepTag]   : fTR->ElPt [indexLepTag],
	  	                  fIsMu ? fTR->MuEta[indexLepTag]   : fTR->ElEta[indexLepTag],
	  	                  fIsMu ? fTR->MuPhi[indexLepTag]   : fTR->ElPhi[indexLepTag],
	  	                  massLep);
	
	lepProbe.SetPtEtaPhiM(fIsMu ? fTR->MuPt [indexLepProbe] : fTR->ElPt [indexLepProbe],
	  		              fIsMu ? fTR->MuEta[indexLepProbe] : fTR->ElEta[indexLepProbe],
	  		              fIsMu ? fTR->MuPhi[indexLepProbe] : fTR->ElPhi[indexLepProbe],
	  		              massLep);
	
	
	fTagPt     = lepTag.Pt();
	fTagEta    = lepTag.Eta();
	fTagPhi    = lepTag.Phi();
	fTagIsoRel = fIsMu ? MuPFIso     (indexLepTag) : ElPFIso     (indexLepTag);
	fTagD0     = fIsMu ? fTR->MuD0PV [indexLepTag] : fTR->ElD0PV [indexLepTag];
	fTagDz     = fIsMu ? fTR->MuDzPV [indexLepTag] : fTR->ElDzPV [indexLepTag];
	fTagPassID = fIsMu ? PassOnlyIDMu(indexLepTag) : PassOnlyIDEl(indexLepTag);

	fTisMuEvent = fIsMu;
	
	fProbePt     = lepProbe.Pt();
	fProbeEta    = lepProbe.Eta();
	fProbePhi    = lepProbe.Phi();
	fProbeIsoRel = fIsMu ? MuPFIso     (indexLepProbe) : ElPFIso     (indexLepProbe);
	fProbeD0     = fIsMu ? fTR->MuD0PV [indexLepProbe] : fTR->ElD0PV [indexLepProbe];
	fProbeDz     = fIsMu ? fTR->MuDzPV [indexLepProbe] : fTR->ElDzPV [indexLepProbe];
	fProbePassID = fIsMu ? PassOnlyIDMu(indexLepProbe) : PassOnlyIDEl(indexLepProbe);
	
	
	fMass2L  =  (lepTag+lepProbe).M();
	fDeltaR  =  lepTag.DeltaR(lepProbe);
	fNvtx    =  fTR->NVrtx;
	fRho     =  fTR->Rho;
	fPfMet   =  fTR->PFType1MET;
	fHt      =  QuickHT();
	if(!fIsData)
		fNTrueInt = fTR->PUnumInteractions;
	else
		fNTrueInt = -99;
	
	// ----
	fAnalysisTree->Fill();
}

bool TnPAnalysis::PassOnlyIDMu(int index){

	if(fTR->MuIsPFMuon[index] == 0)      return false;
	if(fTR->MuIsGlobalMuon[index] == 0)  return false;
	if(fTR->MuNChi2[index] > 10)         return false;
	if(fTR->MuNSiLayers[index] < 6)      return false;
	if(fTR->MuNPxHits[index] < 1)        return false;
	if(fTR->MuNGlMuHits [index] < 1)     return false;
	if(fTR->MuNMatchedStations[index] < 2)   return false;
	if(fTR->MuIso03EMVetoEt[index] > 4.0)  return false;
	if(fTR->MuIso03HadVetoEt[index] > 6.0) return false;
	return true;

}

bool TnPAnalysis::PassOnlyIDEl(int index){

	// this is the tight ID working point:
	if( fabs(fTR->ElSCEta[index]) < 1.479 ){ // Barrel
	    if(fTR->ElSigmaIetaIeta                 [index]  > 0.01 ) return false;
	    if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index]) > 0.06 ) return false;
	    if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index]) > 0.004) return false;
	    if(fTR->ElHcalOverEcal                  [index]  > 0.10 ) return false; // tightened to trigger cut
	}
	if( fabs(fTR->ElSCEta[index]) > 1.479 ){ // Endcap
	    if(fTR->ElSigmaIetaIeta                 [index]  > 0.03 ) return false;
	    if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index]) > 0.03 ) return false;
	    if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index]) > 0.007) return false;
	    if(fTR->ElHcalOverEcal                  [index]  > 0.075 )return false;
	}
	if(fabs(1/fTR->ElCaloEnergy[index] - fTR->ElESuperClusterOverP[index]/fTR->ElCaloEnergy[index]) > 0.05 ) return false;
	if(fTR->ElNumberOfMissingInnerHits[index] > 0 ) return false;
	if(!fTR->ElPassConversionVeto[index] )          return false;
	if(fTR->ElCInfoIsGsfCtfScPixCons[index] == 0 )  return false;
	
	// GAP VETO
	if (fabs(fTR->ElSCEta[index]) > 1.4442 && fabs(fTR->ElSCEta[index]) < 1.566) return false;

	return true;
}


float TnPAnalysis::QuickHT(){
	float ht(0.);

	vector<int> selectedJetInd = PFJetSelection(15., 2.5, &UserAnalysisBase::IsGoodBasicPFJet);
	int njets = std::min( (int)selectedJetInd.size(), 20);

	for(int ind = 0; ind < njets; ind++){
		int jetindex = selectedJetInd[ind];
		float pt = fTR->JPt[jetindex];
		if(pt > 30.) ht += pt;
	}

	return ht;
}
