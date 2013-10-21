#include "DiPhotonJetsAnalyzer.hh"
#include "DiPhotonPurity.hh"
#include "DiPhotonMiniTree.hh"


#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

DiPhotonJetsAnalyzer::DiPhotonJetsAnalyzer(TTree *tree, std::string dataType, Float_t aw, Float_t* _kfac, Float_t _minthrpfphotoncandEB, Float_t _minthrpfphotoncandEE, bool _isstep2, TString _input_filename, UInt_t _uuid) : TreeAnalyzerBase(tree), AddWeight(aw), kfactors(_kfac), minthrpfphotoncandEB(_minthrpfphotoncandEB), minthrpfphotoncandEE(_minthrpfphotoncandEE), isstep2(_isstep2), input_filename(_input_filename), uuid(_uuid) {
  //  fDiPhotonPurity = new DiPhotonPurity(fTR,dataType,AddWeight);

  /*
fTR->fChain->SetBranchStatus("*",0);
fTR->fChain->SetBranchStatus("Run",1);
fTR->fChain->SetBranchStatus("Event",1);
fTR->fChain->SetBranchStatus("LumiSection",1);
fTR->fChain->SetBranchStatus("PUnumInteractions",1);
fTR->fChain->SetBranchStatus("PUnumTrueInteractions",1);
fTR->fChain->SetBranchStatus("PUOOTnumInteractionsEarly",1);
fTR->fChain->SetBranchStatus("PUOOTnumInteractionsLate",1);
fTR->fChain->SetBranchStatus("Rho",1);
fTR->fChain->SetBranchStatus("Sigma",1);
fTR->fChain->SetBranchStatus("Weight",1);
fTR->fChain->SetBranchStatus("HLTResults",1);
fTR->fChain->SetBranchStatus("HLTPrescale",1);
fTR->fChain->SetBranchStatus("PrimVtxGood",1);
fTR->fChain->SetBranchStatus("PrimVtxx",1);
fTR->fChain->SetBranchStatus("PrimVtxy",1);
fTR->fChain->SetBranchStatus("PrimVtxz",1);
fTR->fChain->SetBranchStatus("PrimVtxNChi2",1);
fTR->fChain->SetBranchStatus("PrimVtxNdof",1);
fTR->fChain->SetBranchStatus("PrimVtxIsFake",1);
fTR->fChain->SetBranchStatus("PrimVtxPtSum",1);
fTR->fChain->SetBranchStatus("GoodEvent",1);
fTR->fChain->SetBranchStatus("HBHENoiseFlag",1);
fTR->fChain->SetBranchStatus("HBHENoiseFlagIso",1);
fTR->fChain->SetBranchStatus("CSCTightHaloID",1);
fTR->fChain->SetBranchStatus("EcalDeadTPFilterFlag",1);
fTR->fChain->SetBranchStatus("RA2TrackingFailureFilterFlag",1);
fTR->fChain->SetBranchStatus("NGenPhotons",1);
fTR->fChain->SetBranchStatus("GenPhotonPt",1);
fTR->fChain->SetBranchStatus("GenPhotonEta",1);
fTR->fChain->SetBranchStatus("GenPhotonPhi",1);
fTR->fChain->SetBranchStatus("GenPhotonMotherID",1);
fTR->fChain->SetBranchStatus("GenPhotonMotherStatus",1);
fTR->fChain->SetBranchStatus("GenPhotonIsoDR04",1);
fTR->fChain->SetBranchStatus("NVrtx",1);
fTR->fChain->SetBranchStatus("VrtxX",1);
fTR->fChain->SetBranchStatus("VrtxY",1);
fTR->fChain->SetBranchStatus("VrtxZ",1);
fTR->fChain->SetBranchStatus("VrtxNdof",1);
fTR->fChain->SetBranchStatus("NMus",1);
fTR->fChain->SetBranchStatus("NMusTot",1);
fTR->fChain->SetBranchStatus("MuGood",1);
fTR->fChain->SetBranchStatus("MuIsIso",1);
fTR->fChain->SetBranchStatus("MuIsGlobalMuon",1);
fTR->fChain->SetBranchStatus("MuIsTrackerMuon",1);
fTR->fChain->SetBranchStatus("MuPx",1);
fTR->fChain->SetBranchStatus("MuPy",1);
fTR->fChain->SetBranchStatus("MuPz",1);
fTR->fChain->SetBranchStatus("MuPt",1);
fTR->fChain->SetBranchStatus("MuPtE",1);
fTR->fChain->SetBranchStatus("MuE",1);
fTR->fChain->SetBranchStatus("MuEt",1);
fTR->fChain->SetBranchStatus("MuEta",1);
fTR->fChain->SetBranchStatus("MuPhi",1);
fTR->fChain->SetBranchStatus("MuCharge",1);
fTR->fChain->SetBranchStatus("MuRelIso03",1);
fTR->fChain->SetBranchStatus("MuIso03SumPt",1);
fTR->fChain->SetBranchStatus("MuIso03EmEt",1);
fTR->fChain->SetBranchStatus("MuIso03HadEt",1);
fTR->fChain->SetBranchStatus("MuIso03EMVetoEt",1);
fTR->fChain->SetBranchStatus("MuIso03HadVetoEt",1);
fTR->fChain->SetBranchStatus("MuEem",1);
fTR->fChain->SetBranchStatus("MuEhad",1);
fTR->fChain->SetBranchStatus("MuD0PV",1);
fTR->fChain->SetBranchStatus("MuD0E",1);
fTR->fChain->SetBranchStatus("MuDzPV",1);
fTR->fChain->SetBranchStatus("MuDzE",1);
fTR->fChain->SetBranchStatus("MuNChi2",1);
fTR->fChain->SetBranchStatus("MuNGlHits",1);
fTR->fChain->SetBranchStatus("MuNMuHits",1);
fTR->fChain->SetBranchStatus("MuNTkHits",1);
fTR->fChain->SetBranchStatus("MuNPxHits",1);
fTR->fChain->SetBranchStatus("MuNMatches",1);
fTR->fChain->SetBranchStatus("NEles",1);
fTR->fChain->SetBranchStatus("ElGood",1);
fTR->fChain->SetBranchStatus("ElIsIso",1);
fTR->fChain->SetBranchStatus("ElPx",1);
fTR->fChain->SetBranchStatus("ElPy",1);
fTR->fChain->SetBranchStatus("ElPz",1);
fTR->fChain->SetBranchStatus("ElPt",1);
fTR->fChain->SetBranchStatus("ElE",1);
fTR->fChain->SetBranchStatus("ElEt",1);
fTR->fChain->SetBranchStatus("ElEta",1);
fTR->fChain->SetBranchStatus("ElSCEta",1);
fTR->fChain->SetBranchStatus("ElPhi",1);
fTR->fChain->SetBranchStatus("ElD0PV",1);
fTR->fChain->SetBranchStatus("ElD0E",1);
fTR->fChain->SetBranchStatus("ElDzPV",1);
fTR->fChain->SetBranchStatus("ElDzE",1);
fTR->fChain->SetBranchStatus("ElRelIso03",1);
fTR->fChain->SetBranchStatus("ElDR03TkSumPt",1);
fTR->fChain->SetBranchStatus("ElDR03EcalRecHitSumEt",1);
fTR->fChain->SetBranchStatus("ElDR03HcalTowerSumEt",1);
fTR->fChain->SetBranchStatus("ElCharge",1);
fTR->fChain->SetBranchStatus("ElCInfoIsGsfCtfCons",1);
fTR->fChain->SetBranchStatus("ElScPixCharge",1);
fTR->fChain->SetBranchStatus("ElIDsimpleWP80relIso",1);
fTR->fChain->SetBranchStatus("ElIDsimpleWP95relIso",1);
fTR->fChain->SetBranchStatus("ElEcalDriven",1);
fTR->fChain->SetBranchStatus("ElTrackerDriven",1);
fTR->fChain->SetBranchStatus("Elfbrem",1);
fTR->fChain->SetBranchStatus("ElHcalOverEcal",1);
fTR->fChain->SetBranchStatus("ElSigmaIetaIeta",1);
fTR->fChain->SetBranchStatus("ElDeltaPhiSuperClusterAtVtx",1);
fTR->fChain->SetBranchStatus("ElDeltaEtaSuperClusterAtVtx",1);
fTR->fChain->SetBranchStatus("ElESuperClusterOverP",1);
fTR->fChain->SetBranchStatus("ElNumberOfMissingInnerHits",1);
fTR->fChain->SetBranchStatus("ElConvPartnerTrkDist",1);
fTR->fChain->SetBranchStatus("ElConvPartnerTrkDCot",1);
fTR->fChain->SetBranchStatus("ElConvPartnerTrkCharge",1);
fTR->fChain->SetBranchStatus("NPfCand",1);
fTR->fChain->SetBranchStatus("PfCandPdgId",1);
fTR->fChain->SetBranchStatus("PfCandEta",1);
fTR->fChain->SetBranchStatus("PfCandPhi",1);
fTR->fChain->SetBranchStatus("PfCandPx",1);
fTR->fChain->SetBranchStatus("PfCandPy",1);
fTR->fChain->SetBranchStatus("PfCandPz",1);
fTR->fChain->SetBranchStatus("PfCandEnergy",1);
fTR->fChain->SetBranchStatus("PfCandPt",1);
fTR->fChain->SetBranchStatus("PfCandVx",1);
fTR->fChain->SetBranchStatus("PfCandVy",1);
fTR->fChain->SetBranchStatus("PfCandVz",1);
fTR->fChain->SetBranchStatus("PfCandHasHitInFirstPixelLayer",1);
fTR->fChain->SetBranchStatus("PfCandTrackRefPx",1);
fTR->fChain->SetBranchStatus("PfCandTrackRefPy",1);
fTR->fChain->SetBranchStatus("PfCandTrackRefPz",1);
fTR->fChain->SetBranchStatus("NPhotons",1);
fTR->fChain->SetBranchStatus("PhoGood",1);
fTR->fChain->SetBranchStatus("PhoIsIso",1);
fTR->fChain->SetBranchStatus("PhoPt",1);
fTR->fChain->SetBranchStatus("PhoPx",1);
fTR->fChain->SetBranchStatus("PhoPy",1);
fTR->fChain->SetBranchStatus("PhoPz",1);
fTR->fChain->SetBranchStatus("PhoEta",1);
fTR->fChain->SetBranchStatus("PhoPhi",1);
fTR->fChain->SetBranchStatus("PhoEnergy",1);
fTR->fChain->SetBranchStatus("PhoIso03Ecal",1);
fTR->fChain->SetBranchStatus("PhoIso03Hcal",1);
fTR->fChain->SetBranchStatus("PhoIso03TrkSolid",1);
fTR->fChain->SetBranchStatus("PhoIso03TrkHollow",1);
fTR->fChain->SetBranchStatus("PhoIso03",1);
fTR->fChain->SetBranchStatus("PhoIso04Ecal",1);
fTR->fChain->SetBranchStatus("PhoIso04Hcal",1);
fTR->fChain->SetBranchStatus("PhoIso04TrkSolid",1);
fTR->fChain->SetBranchStatus("PhoIso04TrkHollow",1);
fTR->fChain->SetBranchStatus("PhoIso04",1);
fTR->fChain->SetBranchStatus("PhoHoverE",1);
fTR->fChain->SetBranchStatus("PhoSigmaIetaIeta",1);
fTR->fChain->SetBranchStatus("PhoHasPixSeed",1);
fTR->fChain->SetBranchStatus("PhoPassConvSafeElectronVeto",1);
fTR->fChain->SetBranchStatus("PhoHasConvTrks",1);
fTR->fChain->SetBranchStatus("PhoScSeedSeverity",1);
fTR->fChain->SetBranchStatus("PhoE1OverE9",1);
fTR->fChain->SetBranchStatus("PhoS4OverS1",1);
fTR->fChain->SetBranchStatus("PhoSigmaEtaEta",1);
fTR->fChain->SetBranchStatus("PhoE1x5",1);
fTR->fChain->SetBranchStatus("PhoE2x5",1);
fTR->fChain->SetBranchStatus("PhoE3x3",1);
fTR->fChain->SetBranchStatus("PhoE5x5",1);
fTR->fChain->SetBranchStatus("PhomaxEnergyXtal",1);
fTR->fChain->SetBranchStatus("PhoIso03HcalDepth1",1);
fTR->fChain->SetBranchStatus("PhoIso03HcalDepth2",1);
fTR->fChain->SetBranchStatus("PhoIso04HcalDepth1",1);
fTR->fChain->SetBranchStatus("PhoIso04HcalDepth2",1);
fTR->fChain->SetBranchStatus("PhoIso03nTrksSolid",1);
fTR->fChain->SetBranchStatus("PhoIso03nTrksHollow",1);
fTR->fChain->SetBranchStatus("PhoIso04nTrksSolid",1);
fTR->fChain->SetBranchStatus("PhoIso04nTrksHollow",1);
fTR->fChain->SetBranchStatus("PhoisEB",1);
fTR->fChain->SetBranchStatus("PhoisEE",1);
fTR->fChain->SetBranchStatus("PhoMCmatchindex",1);
fTR->fChain->SetBranchStatus("PhoMCmatchexitcode",1);
fTR->fChain->SetBranchStatus("Pho_ChargedHadronIso",1);
fTR->fChain->SetBranchStatus("Pho_NeutralHadronIso",1);
fTR->fChain->SetBranchStatus("Pho_PhotonIso",1);
fTR->fChain->SetBranchStatus("Pho_isPFPhoton",1);
fTR->fChain->SetBranchStatus("Pho_isPFElectron",1);
fTR->fChain->SetBranchStatus("PhotSCindex",1);
fTR->fChain->SetBranchStatus("pho_matchedPFPhotonCand",1);
fTR->fChain->SetBranchStatus("PhoVx",1);
fTR->fChain->SetBranchStatus("PhoVy",1);
fTR->fChain->SetBranchStatus("PhoVz",1);
fTR->fChain->SetBranchStatus("pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx",1);
fTR->fChain->SetBranchStatus("pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx",1);
fTR->fChain->SetBranchStatus("pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx",1);
fTR->fChain->SetBranchStatus("pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx",1);
fTR->fChain->SetBranchStatus("pho_Cone01NeutralHadronIso_mvVtx",1);
fTR->fChain->SetBranchStatus("pho_Cone02NeutralHadronIso_mvVtx",1);
fTR->fChain->SetBranchStatus("pho_Cone03NeutralHadronIso_mvVtx",1);
fTR->fChain->SetBranchStatus("pho_Cone04NeutralHadronIso_mvVtx",1);
fTR->fChain->SetBranchStatus("pho_Cone01ChargedHadronIso_dR02_dz02_dxy01",1);
fTR->fChain->SetBranchStatus("pho_Cone02ChargedHadronIso_dR02_dz02_dxy01",1);
fTR->fChain->SetBranchStatus("pho_Cone03ChargedHadronIso_dR02_dz02_dxy01",1);
fTR->fChain->SetBranchStatus("pho_Cone04ChargedHadronIso_dR02_dz02_dxy01",1);
fTR->fChain->SetBranchStatus("pho_Cone03PFCombinedIso",1);
fTR->fChain->SetBranchStatus("pho_Cone04PFCombinedIso",1);
fTR->fChain->SetBranchStatus("NSuperClusters",1);
fTR->fChain->SetBranchStatus("SCRaw",1);
fTR->fChain->SetBranchStatus("SCPre",1);
fTR->fChain->SetBranchStatus("SCEta",1);
fTR->fChain->SetBranchStatus("SCPhi",1);
fTR->fChain->SetBranchStatus("SCPhiWidth",1);
fTR->fChain->SetBranchStatus("SCEtaWidth",1);
fTR->fChain->SetBranchStatus("SCBrem",1);
fTR->fChain->SetBranchStatus("SCR9",1);
fTR->fChain->SetBranchStatus("SCx",1);
fTR->fChain->SetBranchStatus("SCy",1);
fTR->fChain->SetBranchStatus("SCz",1);
fTR->fChain->SetBranchStatus("SCNXtals",1);
fTR->fChain->SetBranchStatus("SCxtalX",1);
fTR->fChain->SetBranchStatus("SCxtalY",1);
fTR->fChain->SetBranchStatus("SCxtalZ",1);
fTR->fChain->SetBranchStatus("SCxtalEtaWidth",1);
fTR->fChain->SetBranchStatus("SCxtalPhiWidth",1);
fTR->fChain->SetBranchStatus("SCxtalfrontX",1);
fTR->fChain->SetBranchStatus("SCxtalfrontY",1);
fTR->fChain->SetBranchStatus("SCxtalfrontZ",1);
fTR->fChain->SetBranchStatus("NJets",1);
fTR->fChain->SetBranchStatus("JGood",1);
fTR->fChain->SetBranchStatus("JPx",1);
fTR->fChain->SetBranchStatus("JPy",1);
fTR->fChain->SetBranchStatus("JPz",1);
fTR->fChain->SetBranchStatus("JPt",1);
fTR->fChain->SetBranchStatus("JE",1);
fTR->fChain->SetBranchStatus("JEt",1);
fTR->fChain->SetBranchStatus("JEta",1);
fTR->fChain->SetBranchStatus("JPhi",1);
fTR->fChain->SetBranchStatus("JEcorr",1);
//fTR->fChain->SetBranchStatus("JArea",1);
//fTR->fChain->SetBranchStatus("JNConstituents",1);
//fTR->fChain->SetBranchStatus("JNAssoTracks",1);
//fTR->fChain->SetBranchStatus("JChargedEmFrac",1);
//fTR->fChain->SetBranchStatus("JNeutralEmFrac",1);
//fTR->fChain->SetBranchStatus("JChargedHadFrac",1);
//fTR->fChain->SetBranchStatus("JNeutralHadFrac",1);
//fTR->fChain->SetBranchStatus("JbTagProbTkCntHighEff",1);
//fTR->fChain->SetBranchStatus("JbTagProbTkCntHighPur",1);
//fTR->fChain->SetBranchStatus("JbTagProbSimpSVHighEff",1);
//fTR->fChain->SetBranchStatus("JbTagProbSimpSVHighPur",1);
//fTR->fChain->SetBranchStatus("CAJE",1);
//fTR->fChain->SetBranchStatus("CAJEt",1);
//fTR->fChain->SetBranchStatus("CAJEta",1);
//fTR->fChain->SetBranchStatus("CAJID_HPD",1);
//fTR->fChain->SetBranchStatus("CAJID_RBX",1);
//fTR->fChain->SetBranchStatus("CAJID_n90Hits",1);
//fTR->fChain->SetBranchStatus("CAJEMfrac",1);
//fTR->fChain->SetBranchStatus("CAJChfrac",1);
//fTR->fChain->SetBranchStatus("PF2PAT3NJets",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JPx",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JPy",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JPz",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JPt",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JE",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JEt",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JEta",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JPhi",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JbTagProbTkCntHighEff",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JbTagProbTkCntHighPur",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JbTagProbSimpSVHighEff",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JbTagProbSimpSVHighPur",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JIDLoose",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JChMult",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JNeuMult",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JChHadfrac",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JNeuHadfrac",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JChEmfrac",1);
//fTR->fChain->SetBranchStatus("PF2PAT3JNeuEmfrac",1);
//fTR->fChain->SetBranchStatus("PfTau3NObjs",1);
//fTR->fChain->SetBranchStatus("PfTau3Px",1);
//fTR->fChain->SetBranchStatus("PfTau3Py",1);
//fTR->fChain->SetBranchStatus("PfTau3Pz",1);
//fTR->fChain->SetBranchStatus("PfTau3E",1);
//fTR->fChain->SetBranchStatus("PfTau3DecayMode",1);
//fTR->fChain->SetBranchStatus("SumEt",1);
//fTR->fChain->SetBranchStatus("RawMET",1);
//fTR->fChain->SetBranchStatus("RawMETSignificance",1);
//fTR->fChain->SetBranchStatus("TCMET",1);
//fTR->fChain->SetBranchStatus("TCMETphi",1);
//fTR->fChain->SetBranchStatus("MuJESCorrMET",1);
//fTR->fChain->SetBranchStatus("MuJESCorrMETphi",1);
//fTR->fChain->SetBranchStatus("PFMET",1);
//fTR->fChain->SetBranchStatus("PFMETphi",1);
*/

  fDiPhotonMiniTree = new DiPhotonMiniTree(fTR,dataType,AddWeight,kfactors,minthrpfphotoncandEB,minthrpfphotoncandEE,isstep2,input_filename,uuid);
}

DiPhotonJetsAnalyzer::~DiPhotonJetsAnalyzer(){
  //	delete fDiPhotonPurity;
	delete fDiPhotonMiniTree;
	if(!fTR->fChain) cout << "DiPhotonJetsAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void DiPhotonJetsAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;

	if(fMaxEvents==-1)nentries=fTR->GetEntries();
	if(fMaxEvents>0)nentries=fMaxEvents;
	
	// loop over all ntuple entries
	// nentries = 200;
	for( Long64_t jentry = 0; jentry < nentries; jentry++ ){
		PrintProgress(jentry);
		fTR->GetEntry(jentry);
		if ( fCurRun != fTR->Run ) {
		  fCurRun = fTR->Run;
		  //		  fDiPhotonPurity->BeginRun(fCurRun);
		  fDiPhotonMiniTree->BeginRun(fCurRun);
		  skipRun = false;
		  if ( !CheckRun() ) skipRun = true;
		}
		// Check if new lumi is in JSON file
		if ( !skipRun && fCurLumi != fTR->LumiSection ) {
		  fCurLumi = fTR->LumiSection;
		  skipLumi = false; // Re-initialise
		  if ( !CheckRunLumi() ) skipLumi = true;
		}
		if ( !(skipRun || skipLumi) ) {
		  //		  fDiPhotonPurity->Analyze();
		  fDiPhotonMiniTree->Analyze();
		}
	}
}

// Method called before starting the event loop
void DiPhotonJetsAnalyzer::BeginJob(string fdata_PileUp, string fmc_PileUp){
// 	fDiPhotonPurity->SetOutputDir(fOutputDir);
// 	fDiPhotonPurity->SetOutputFile(fOutputFile);
// 	fDiPhotonPurity->fVerbose = fVerbose;
// 	fDiPhotonPurity->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
// 	fDiPhotonPurity->Begin();

	fDiPhotonMiniTree->SetOutputDir(fOutputDir);
	fDiPhotonMiniTree->SetOutputFile(fOutputFile);
	fDiPhotonMiniTree->fVerbose = fVerbose;
	fDiPhotonMiniTree->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
	fDiPhotonMiniTree->SetPileUp3DSrc(fdata_PileUp, fmc_PileUp);
	fDiPhotonMiniTree->Begin();

}

// Method called after finishing the event loop
void DiPhotonJetsAnalyzer::EndJob(){
  //	fDiPhotonPurity->End();
	fDiPhotonMiniTree->End();
}

