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

	ReadPDGTable("/shome/mdunser/ttW2013/CMSSW_5_3_7_patch5/src/ASAnalysis/pdgtable.txt");
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
	fAnalysisTree->Branch("tag3DIP",            &fTag3DIP,         "tag3DIP/F");
	fAnalysisTree->Branch("tagPassID",        &fTagPassID,     "tagPassID/F");
	//
	fAnalysisTree->Branch("probePt",            &fProbePt,         "probePt/F");
	fAnalysisTree->Branch("probeEta",           &fProbeEta,        "probeEta/F");
	fAnalysisTree->Branch("probePhi",           &fProbePhi,        "probePhi/F");
	fAnalysisTree->Branch("probePtErr",         &fProbePtErr,      "probePtErr/F");
	fAnalysisTree->Branch("probeIsoRel",        &fProbeIsoRel,     "probeIsoRel/F");
	fAnalysisTree->Branch("probeD0",            &fProbeD0,         "probeD0/F");
	fAnalysisTree->Branch("probeDz",            &fProbeDz,         "probeDz/F");
	fAnalysisTree->Branch("probe3DIP",            &fProbe3DIP,         "probe3DIP/F");
	fAnalysisTree->Branch("probePassID",        &fProbePassID,     "probePassID/F");
	fAnalysisTree->Branch("probeIsGlobal",      &fProbeIsGlobal,   "probeIsGlobal/I");
	
	fAnalysisTree->Branch("mass2L",        &fMass2L,     "mass2L/F");
	fAnalysisTree->Branch("mllMatchedLep", &fMass2LMatchedLeptons,  "mllMatchedLep/F");
	fAnalysisTree->Branch("mGenZ",         &fMass2LFromZ,     "mGenZ/F");
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
	fTisMuEvent = fIsMu;
	
	
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
	
	
	// TAG VARIABLES
	fTagPt     = lepTag.Pt();
	fTagEta    = lepTag.Eta();
	fTagPhi    = lepTag.Phi();
	fTagIsoRel = fIsMu ? MuPFIso       (indexLepTag) : ElPFIso       (indexLepTag);
	fTagD0     = fIsMu ? fTR->MuD0PV   [indexLepTag] : fTR->ElD0PV   [indexLepTag];
	fTagDz     = fIsMu ? fTR->MuDzPV   [indexLepTag] : fTR->ElDzPV   [indexLepTag];
	fTag3DIP   = fIsMu ? fTR->MuD03DPV [indexLepTag] : fTR->ElD03DPV [indexLepTag];
	fTagPassID = fIsMu ? PassOnlyIDMu  (indexLepTag) : PassOnlyIDEl  (indexLepTag);

	
	// PROBE VARIABLES
	fProbePt     = lepProbe.Pt();
	fProbeEta    = lepProbe.Eta();
	fProbePhi    = lepProbe.Phi();
	fProbePtErr  = fIsMu ? fTR->MuPtE    [indexLepProbe] : fTR->ElPtE    [indexLepProbe];
	fProbeIsoRel = fIsMu ? MuPFIso       (indexLepProbe) : ElPFIso       (indexLepProbe);
	fProbeD0     = fIsMu ? fTR->MuD0PV   [indexLepProbe] : fTR->ElD0PV   [indexLepProbe];
	fProbeDz     = fIsMu ? fTR->MuDzPV   [indexLepProbe] : fTR->ElDzPV   [indexLepProbe];
	fProbe3DIP   = fIsMu ? fTR->MuD03DPV [indexLepProbe] : fTR->ElD03DPV [indexLepProbe];
	fProbePassID = fIsMu ? PassOnlyIDMu  (indexLepProbe) : PassOnlyIDEl  (indexLepProbe);
	fProbeIsGlobal = fIsMu ? fTR->MuIsGlobalMuon[indexLepProbe] : -1;
	
	
	fMass2L  =  (lepTag+lepProbe).M();
	// get the generator mass if the leptons are matched to a signal lepton
	if (!fIsData) {
		float tPt, tEta, tPhi;
		float pPt, pEta, pPhi;
		bool tagMatch   = fIsMu ? IsSignalMuon(indexLepTag, tPt, tEta, tPhi)   : IsSignalElectron(indexLepTag, tPt, tEta, tPhi) ;
		bool probeMatch = fIsMu ? IsSignalMuon(indexLepProbe, pPt, pEta, pPhi) : IsSignalElectron(indexLepProbe, pPt, pEta, pPhi) ;
		TLorentzVector lep1, lep2;
		
		if (tagMatch && probeMatch){
			lep1.SetPtEtaPhiM(tPt, tEta, tPhi, massLep );
			lep2.SetPtEtaPhiM(pPt, pEta, pPhi, massLep );
			fMass2LMatchedLeptons = (lep1+lep2).M();
		// cout << " mass: " << (lep1+lep2).M() << endl;
		}
		else fMass2LMatchedLeptons = -1.;

		// get the mass from the status 3 Z:
		fMass2LFromZ = getGenZMass();
	}
	else{
		fMass2LMatchedLeptons = -99.;
		fMass2LFromZ          = -99.;
	}


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


// ----------------------------
// --  MATCHING STUFF  --------
// ----------------------------

bool TnPAnalysis::IsSignalElectron(int index, float &elpt, float &eleta, float &elphi){

	int matched = -1;
	// Match to a gen electron
	float mindr(100.);
	for(size_t i = 0; i < fTR->NGenLeptons; ++i){
		if(abs(fTR->GenLeptonID[i]) != 11) continue; // electrons
		float eta = fTR->GenLeptonEta[i];
		float phi = fTR->GenLeptonPhi[i];
		float DR = Util::GetDeltaR(eta, fTR->ElEta[index], phi, fTR->ElPhi[index]);
		if(DR > mindr) continue; // minimize DR
		mindr = DR;
		matched = i;
	}

	bool isGenElectron = true;
        bool isSignalGenElectron = false;	
	if(mindr > 0.2) isGenElectron = false; // max distance
	if(matched < 0) isGenElectron = false; // match unsuccessful

	pdgparticle mo, gmo;
	if (isGenElectron){
	  GetPDGParticle(mo,  abs(fTR->GenLeptonMID [matched]));
	  GetPDGParticle(gmo, abs(fTR->GenLeptonGMID[matched]));
	  if(mo.get_type() > 10  || gmo.get_type() > 10)  isSignalGenElectron = false; // mother is SM hadron
	  if(mo.get_type() == 9  || gmo.get_type() == 9)  isSignalGenElectron = true; // matched to susy
	  if(mo.get_type() == 4 && abs(mo.get_pdgid()) != 21 && abs(mo.get_pdgid()) != 22) isSignalGenElectron = true; // mother is gauge or higgs, but not gamma or gluon
	  if(abs(mo.get_pdgid()) == 15 && (abs(gmo.get_pdgid()) == 6 ||abs(gmo.get_pdgid()) == 23 ||abs(gmo.get_pdgid()) == 24 ) ) isSignalGenElectron = true; // mother is a tau and grandma is a top/W/Z
	  if(abs(mo.get_pdgid()) == 6) isSignalGenElectron = true; // mother is a top
	  // If we've matched to a gen electron, then we'll set the ids now and exit
	  elpt  = fTR->GenLeptonPt [matched];
	  eleta = fTR->GenLeptonEta[matched];
	  elphi = fTR->GenLeptonPhi[matched];
	  return isSignalGenElectron;
	}

	// Otherwise, we'll look for a gen particle in a larger(?) cone
	// and match to whatever is closest
	int matchedPart = -1;
	// Match to a gen particle
	float mindrPart(100.);
	for(size_t i = 0; i < fTR->nGenParticles; ++i){
	  if(fTR->genInfoStatus[i] != 1  && fTR->genInfoStatus[i] != 3 ) continue; 
	  float eta = fTR->genInfoEta[i];
	  float phi = fTR->genInfoPhi[i];
	  float DR = Util::GetDeltaR(eta, fTR->ElEta[index], phi, fTR->ElPhi[index]);
	  if(DR > mindrPart) continue; // minimize DR
	  mindrPart = DR;
	  matchedPart = i;
	}

        bool isGenPart = true;
        bool isSignalGenPart = false;
	if(mindrPart > 0.3) isGenPart = false; // max distance
	if(matchedPart < 0) isGenPart = false; // match unsuccessful

	if (isGenPart){
	  //keep going back until we get a mother with a different id
	  int moIndex = fTR->genInfoMo1[matchedPart];
	  while (fTR->genInfoId[matchedPart] == fTR->genInfoId[moIndex]){
	    moIndex  = fTR->genInfoMo1[moIndex];
	  }
	  //keep going back until we get a grandmother with a different id
	  int gmoIndex = fTR->genInfoMo1[moIndex];
	  while (fTR->genInfoId[moIndex] == fTR->genInfoId[gmoIndex]){
	    gmoIndex  = fTR->genInfoMo1[gmoIndex];
	  }
	  GetPDGParticle(mo ,  abs(fTR->genInfoId[moIndex] ));
	  GetPDGParticle(gmo,  abs(fTR->genInfoId[gmoIndex]));	
	  elpt  = fTR->genInfoPt [matchedPart];
	  eleta = fTR->genInfoEta[matchedPart];
	  elphi = fTR->genInfoPhi[matchedPart];
	  return false;
	}

	//If we really didn't find a match to a lepton or any other particle, set everything to zero
	elpt  = 0.;
	eleta = 0.;
	elphi = 0.;
	return false;	
}

bool TnPAnalysis::IsSignalMuon(int index, float &mupt, float &mueta, float &muphi){

	int matched = -1;
	// Match to a gen muon
	float mindr(100.);
	for(size_t i = 0; i < fTR->NGenLeptons; ++i){
		if(abs(fTR->GenLeptonID[i]) != 13) continue; // muons
		float eta = fTR->GenLeptonEta[i];
		float phi = fTR->GenLeptonPhi[i];
		float DR = Util::GetDeltaR(eta, fTR->MuEta[index], phi, fTR->MuPhi[index]);
		if(DR > mindr) continue; // minimize DR
		mindr = DR;
		matched = i;
	}

	bool isGenMuon = true;
	bool isSignalGenMuon = false;	
	if(mindr > 0.2) isGenMuon = false; // max distance
	if(matched < 0) isGenMuon = false; // match unsuccessful

	pdgparticle mo, gmo;
	if (isGenMuon){
	  GetPDGParticle(mo,  abs(fTR->GenLeptonMID [matched]));
	  GetPDGParticle(gmo, abs(fTR->GenLeptonGMID[matched]));
	  if(mo.get_type() > 10  || gmo.get_type() > 10)  isSignalGenMuon = false; // mother is SM hadron
	  if(mo.get_type() == 9  || gmo.get_type() == 9)  isSignalGenMuon = true; // matched to susy
	  if(mo.get_type() == 4 && abs(mo.get_pdgid()) != 21 && abs(mo.get_pdgid()) != 22) isSignalGenMuon = true; // mother is gauge or higgs, but not gamma or gluon
	  if(abs(mo.get_pdgid()) == 15 && (abs(gmo.get_pdgid()) == 6 ||abs(gmo.get_pdgid()) == 23 ||abs(gmo.get_pdgid()) == 24 ) ) isSignalGenMuon = true; // mother is a tau and grandma is a top/W/Z
	  if(abs(mo.get_pdgid()) == 6) isSignalGenMuon = true; // mother is a top
	  // If we've matched to a gen muon, then we'll set the ids now and exit
	  mupt  = fTR->GenLeptonPt [matched];
	  mueta = fTR->GenLeptonEta[matched];
	  muphi = fTR->GenLeptonPhi[matched];
	  return isSignalGenMuon;
	}

	// Otherwise, we'll look for a gen particle in a larger cone and match to whatever is closest
	int matchedPart = -1;
	// Match to a gen particle
	float mindrPart(100.);
	for(size_t i = 0; i < fTR->nGenParticles; ++i){
	  if(fTR->genInfoStatus[i] != 1  && fTR->genInfoStatus[i] != 3 ) continue; 
	  float eta = fTR->genInfoEta[i];
	  float phi = fTR->genInfoPhi[i];
	  float DR = Util::GetDeltaR(eta, fTR->MuEta[index], phi, fTR->MuPhi[index]);
	  if(DR > mindrPart) continue; // minimize DR
	  mindrPart = DR;
	  matchedPart = i;
	}

	bool isGenPart = true;
	bool isSignalGenPart = false;
	if(mindrPart > 0.3) isGenPart = false; // max distance
	if(matchedPart < 0) isGenPart = false; // match unsuccessful

	if (isGenPart){
	  //keep going back until we get a mother with a different id
	  int moIndex = fTR->genInfoMo1[matchedPart];
	  while (fTR->genInfoId[matchedPart] == fTR->genInfoId[moIndex]){
	    moIndex  = fTR->genInfoMo1[moIndex];
	  }
	  //keep going back until we get a grandmother with a different id
	  int gmoIndex = fTR->genInfoMo1[moIndex];
	  while (fTR->genInfoId[moIndex] == fTR->genInfoId[gmoIndex]){
	    gmoIndex  = fTR->genInfoMo1[gmoIndex];
	  }
	  GetPDGParticle(mo ,  abs(fTR->genInfoId[moIndex] ));
	  GetPDGParticle(gmo,  abs(fTR->genInfoId[gmoIndex]));
	  mupt  = fTR->genInfoPt [matchedPart];
	  mueta = fTR->genInfoEta[matchedPart];
	  muphi = fTR->genInfoPhi[matchedPart];
	  return false;
	}

	//If we really didn't find a match to a lepton or any other particle, set everything to zero
	mupt  = 0.;
	mueta = 0.;
	muphi = 0.;
	return false;
}

float TnPAnalysis::getGenZMass(){
	float mass = -1.;
	for(int i=0; i < fTR->nGenParticles; ++i){
		if(abs(fTR->genInfoId[i]) == 23 && fTR->genInfoStatus[i]  == 3){
			mass = fTR->genInfoM[i];
		}
	}
	return mass;
}
