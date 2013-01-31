#include "SSDLAnalysis.hh"
#include "helper/Monitor.hh"
#include "helper/PUWeight.h"

#include "JetCorrectionUncertainty.h"


using namespace std;

const int SSDLAnalysis::fMaxNjets;
const int SSDLAnalysis::fMaxNmus;
const int SSDLAnalysis::fMaxNeles;
const int SSDLAnalysis::nx;
const float SSDLAnalysis::x_values[nx] =  {0.05, 0.5, 0.95};


TString SSDLAnalysis::gBaseDir = "/shome/mdunser/workspace/CMSSW_5_2_5/src/DiLeptonAnalysis/NTupleProducer/macros/";
// TString SSDLAnalysis::gBaseDir = "/shome/stiegerb/Workspace/cmssw/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/";

//____________________________________________________________________________
SSDLAnalysis::SSDLAnalysis(TreeReader *tr): UserAnalysisBase(tr){
	fHLTPaths.clear();

	fCutnames[0] = "All events";
	fCutnames[1] = " ... passes primary Vertex cuts";
	fCutnames[2] = " ... passes triggers (data only)";
	fCutnames[3] = " ... has at least one loose lepton (mc only)";


	fCounter.fill(fCutnames[0], 0.);
	fCounter.fill(fCutnames[1], 0.);
	fCounter.fill(fCutnames[2], 0.);
	fCounter.fill(fCutnames[3], 0.);
}
SSDLAnalysis::~SSDLAnalysis(){
}

//____________________________________________________________________________
void SSDLAnalysis::Begin(const char* filename){
	// cout << "DEBUG: " << __LINE__ << endl;
	ReadTriggers(gBaseDir + "HLTPaths_SSDL_2012.dat");
	ReadPDGTable(gBaseDir + "pdgtable.txt");

	static const int gM0bins(150), gM0min(0), gM0max(3000), gM12bins(50), gM12min(0), gM12max(1000);
	if(!fIsData){
		fMsugraCount = new TH2D("msugra_count", "msugra_count", gM0bins, gM0min+10, gM0max+10, gM12bins, gM12min+10, gM12max+10);
		for (int i=0; i<10; i++) {
			fProcessCount[i] = new TH2D(Form("msugra_count_process%i",i+1), Form("msugra_count_process%i",i+1), gM0bins, gM0min+10, gM0max+10, gM12bins, gM12min+10, gM12max+10);
		}
		// define all x-values for your scan in the header before running this code!
		fRightHandedSlepCountAll    = new TH2D("RightHandedSlepCountAll" , "RightHandedSlepCountAll" , 101 , -5 , 1005 , 101 , -5 , 1005);
		fRightHandedCountAll        = new TH2D("RightHandedCountAll"     , "RightHandedCountAll"     , 101 , -5 , 1005 , 101 , -5 , 1005);
		fTChiSlepSlepCountAll       = new TH2D("TChiSlepSlepCountAll"    , "TChiSlepSlepCountAll"    , 101 , -5 , 1005 , 101 , -5 , 1005);
		fTChiSlepSnuCountAll        = new TH2D("TChiSlepSnuCountAll"     , "TChiSlepSnuCountAll"     , 101 , -5 , 1005 , 101 , -5 , 1005);
		fModelCountAll              = new TH2D("ModelCountAll"           , "ModelCountAll"           , 101 , -5 , 1005 , 101 , -5 , 1005);
		for (int i = 0; i < nx; ++i) {
			fRightHandedSlepCount[i] = new TH2D(Form("RightHandedSlepCount%.0f" , 100*x_values[i]) , Form("RightHandedSlepCount%.0f" , 100*x_values[i]) , 101 , -5 , 1005 , 101 , -5 , 1005);
			fRightHandedCount[i]     = new TH2D(Form("RightHandedCount%.0f"     , 100*x_values[i]) , Form("RightHandedCount%.0f"     , 100*x_values[i]) , 101 , -5 , 1005 , 101 , -5 , 1005);
			fTChiSlepSlepCount[i]    = new TH2D(Form("TChiSlepSlepCount%.0f"    , 100*x_values[i]) , Form("TChiSlepSlepCount%.0f"    , 100*x_values[i]) , 101 , -5 , 1005 , 101 , -5 , 1005);
			fTChiSlepSnuCount[i]     = new TH2D(Form("TChiSlepSnuCount%.0f"     , 100*x_values[i]) , Form("TChiSlepSnuCount%.0f"     , 100*x_values[i]) , 101 , -5 , 1005 , 101 , -5 , 1005);
			fModelCount[i]           = new TH2D(Form("ModelCount%.0f"           , 100*x_values[i]) , Form("ModelCount%.0f"           , 100*x_values[i]) , 101 , -5 , 1005 , 101 , -5 , 1005);
		}
	}
	BookTree();
	fHEvCount = new TH1F("EventCount", "Event Counter", 1, 0., 1.); // count number of generated events
	if(!fIsData && fDoFillEffTree) BookEffTree();
}

//____________________________________________________________________________
void SSDLAnalysis::End(){
	if (!fIsData) {
		fMsugraCount->Write();
		for (int i=0; i<10; i++) {
			fProcessCount[i]->Write();
		}
		fRightHandedSlepCountAll -> Write();
		fRightHandedCountAll     -> Write();
		fTChiSlepSlepCountAll    -> Write();
		fTChiSlepSnuCountAll     -> Write();
		fModelCountAll           -> Write();
		for (int i = 0; i < nx; ++i) {
			fRightHandedSlepCount [i] -> Write();
			fRightHandedCount     [i] -> Write();
			fTChiSlepSlepCount    [i] -> Write();
			fTChiSlepSnuCount     [i] -> Write();
			fModelCount           [i] -> Write();
		}
	}
	fOutputFile->cd();
	fHEvCount->Write();
	fAnalysisTree->Write();
	if(fDoFillEffTree && !fIsData) fLepEffTree->Write();
	fOutputFile->Close();
	fCounter.print();
}

//____________________________________________________________________________
void SSDLAnalysis::ReadTriggers(const char* triggerfile){
	// Read in a bunch of HLT paths to be stored in the mini tree
	ifstream IN(triggerfile);
	char buffer[200];
	char StringValue[100];
	bool ok(false);
	
	
	while( IN.getline(buffer, 200, '\n') ){
		// ok = false;
		if (buffer[0] == '#') continue; // Skip lines commented with '#'
		if( !strcmp(buffer, "PATHSET")){
			HLTPathSet ps;

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Name\t%s", StringValue);
			ps.name = TString(StringValue);
			while( IN.getline(buffer, 200, '\n') && IN.gcount() > 1 && buffer[0] != '#' ){
				// End declaration of paths of a set with either an empty line or a comment
				sscanf(buffer, "Path\t%s", StringValue);
				ps.paths.push_back(string(StringValue));
			}
			fHLTPathSets.push_back(ps);
			ok = true;
		}

		if(!ok) cout << "%% SSDLAnalysis::ReadTriggers ==> ERROR: Reading failed!" << endl;
	}
	if(fVerbose > 0){
		cout << "Adding Trigger path sets:" << endl;
		for(size_t i = 0; i < fHLTPathSets.size(); ++i){
			HLTPathSet ps = fHLTPathSets[i];
			cout << " " << ps.name << endl;
			for(size_t j = 0; j < ps.paths.size(); ++j){
				cout << "   > " << ps.paths[j] << endl;
			}
		}
		cout << "--------------" << endl;
	}

	fHLTResults.resize(fHLTPathSets.size());
	fHLTPrescales.resize(fHLTPathSets.size());
}
void SSDLAnalysis::AddTriggerBranches(){
	for(unsigned int i = 0; i < fHLTPathSets.size(); i++){
		HLTPathSet ps = fHLTPathSets[i];
		TString prescalename = ps.name + "_PS";
		if(AddBranch(ps.name.Data(),      "I", &fHLTResults[i])   == false ) exit(-1);
		if(AddBranch(prescalename.Data(), "I", &fHLTPrescales[i]) == false ) exit(-1);
	}
}
bool SSDLAnalysis::FillTriggers(){
	// Returns OR of trigger results, i.e. true if ANY of them passed
	bool accept = false;
	if(fHLTPathSets.size() == 0 ) return false;

	for(unsigned int i = 0; i < fHLTPathSets.size(); i++){ // loop over path sets
		HLTPathSet ps = fHLTPathSets[i];

		for(size_t j = 0; j < ps.paths.size(); ++j){ // loop over paths
			fHLTResults[i]   = -1;
			fHLTPrescales[i] = -1;
			string path = ps.paths[j];

			if(GetHLTBit(path) == -1) continue; // Bit not found
			bool triggered   = GetHLTResult(path);
			fHLTResults[i]   = triggered ? 1:0;
			fHLTPrescales[i] = GetHLTPrescale(path);
			accept = accept || triggered;
			break; // break loop after the bit is found
		}
	}
	return accept;
}
const bool SSDLAnalysis::AddBranch( const char* name, const char* type, void* address, const char* size ){
	// This is copied from FillerBase.cc
	// Form input
	std::string fullname(name);

	std::string branchType(fullname);
	if ( size ) // Size needs to be pre-fixed
	branchType += "[" + std::string(size) + "]";
	branchType += "/"+std::string(type);

	// Declare branch
	TBranch* b = fAnalysisTree->Branch(fullname.c_str(),address,branchType.c_str());

	return !(b==0); // return 1 if branch was successfully created, 0 otherwise
}

//____________________________________________________________________________
void SSDLAnalysis::BookTree(){
	fOutputFile->cd();
	fAnalysisTree = new TTree("Analysis", "AnalysisTree");

    // run/sample properties
	fAnalysisTree->Branch("Run",              &fTRunNumber,         "Run/I");
	fAnalysisTree->Branch("Event",            &fTEventNumber,       "Event/I");
	fAnalysisTree->Branch("LumiSec",          &fTLumiSection,       "LumiSec/I");

	fAnalysisTree->Branch("m0",            &fTm0,         "m0/F");
	fAnalysisTree->Branch("m12",           &fTm12,        "m12/F");
	fAnalysisTree->Branch("process",       &fTprocess,    "process/I");

	fAnalysisTree->Branch("mGlu",          &fTmGlu,       "mGlu/F");
	fAnalysisTree->Branch("mLSP",          &fTmLSP,       "mLSP/F");
	fAnalysisTree->Branch("isTChiSlepSnu", &fTisTChiSlepSnu,       "isTChiSlepSnu/I");
	fAnalysisTree->Branch("isRightHanded", &fTisRightHanded,       "isRightHanded/I");

	// HLT triggers
	AddTriggerBranches();
	
	// event properties
	fAnalysisTree->Branch("Rho",           &fTrho,       "Rho/F");
	fAnalysisTree->Branch("NVrtx",         &fTnvrtx,     "NVrtx/I");
	fAnalysisTree->Branch("PUWeight",      &fTpuweight,  "PUWeight/F");

	// single-muon properties
	fAnalysisTree->Branch("NMus"          ,&fTnqmus,          "NMus/I");
	fAnalysisTree->Branch("IsSignalMuon"  ,&fTIsSignalMuon,   "IsSignalMuon[NMus]/I");
	fAnalysisTree->Branch("MuPt"          ,&fTmupt,           "MuPt[NMus]/F");
	fAnalysisTree->Branch("MuEta"         ,&fTmueta,          "MuEta[NMus]/F");
	fAnalysisTree->Branch("MuPhi"         ,&fTmuphi,          "MuPhi[NMus]/F");
	fAnalysisTree->Branch("MuCharge"      ,&fTmucharge,       "MuCharge[NMus]/I");
	fAnalysisTree->Branch("MuDetIso"      ,&fTmudetiso,       "MuDetIso[NMus]/F");
	fAnalysisTree->Branch("MuPFIso"       ,&fTmupfiso,        "MuPFIso[NMus]/F");
	fAnalysisTree->Branch("MuPFIso04"     ,&fTmupfiso04,      "MuPFIso04[NMus]/F");
	fAnalysisTree->Branch("MuPFChIso"     ,&fTmupfchiso,      "MuPFChIso[NMus]/F");
	fAnalysisTree->Branch("MuPFNeIso"     ,&fTmupfneiso,      "MuPFNeIso[NMus]/F");
	fAnalysisTree->Branch("MuPFNeIsoUnc"  ,&fTmupfneisounc,   "MuPFNeIsoUnc[NMus]/F");
	fAnalysisTree->Branch("MuRadIso"      ,&fTmuradiso,       "MuRadIso[NMus]/F");
	fAnalysisTree->Branch("MuD0"          ,&fTmud0,           "MuD0[NMus]/F");
	fAnalysisTree->Branch("MuDz"          ,&fTmudz,           "MuDz[NMus]/F");
	fAnalysisTree->Branch("MuEMVetoEt"    ,&fTmuEMVetoEt,     "MuEMVetoEt[NMus]/F");
	fAnalysisTree->Branch("MuHadVetoEt"   ,&fTmuHadVetoEt,    "MuHadVetoEt[NMus]/F");
	fAnalysisTree->Branch("MuPassesTightID",&fTmuPassesTightID,    "MuPassesTightID[NMus]/I");
	fAnalysisTree->Branch("MuPtE"         ,&fTmuptE,          "MuPtE[NMus]/F");
	fAnalysisTree->Branch("MuGenID"       ,&fTmuid,           "MuGenID[NMus]/I");
	fAnalysisTree->Branch("MuGenMID"      ,&fTmumoid,         "MuGenMID[NMus]/I");
	fAnalysisTree->Branch("MuGenGMID"     ,&fTmugmoid,        "MuGenGMID[NMus]/I");
	fAnalysisTree->Branch("MuGenType"     ,&fTmutype,         "MuGenType[NMus]/I");
	fAnalysisTree->Branch("MuGenMType"    ,&fTmumotype,       "MuGenMType[NMus]/I");
	fAnalysisTree->Branch("MuGenGMType"   ,&fTmugmotype,      "MuGenGMType[NMus]/I");
	fAnalysisTree->Branch("MuMT"          ,&fTmuMT,           "MuMT[NMus]/F");

	// single-electron properties
	fAnalysisTree->Branch("NEls",                   &fTnqels,               "NEls/I");
	fAnalysisTree->Branch("IsSignalElectron" ,      &fTIsSignalElectron,    "IsSignalElectron[NEls]/I");
	fAnalysisTree->Branch("ElCharge",               &fTElcharge,            "ElCh[NEls]/I");
	fAnalysisTree->Branch("ElChIsCons",             &fTElChargeIsCons,      "ElChIsCons[NEls]/I");
	fAnalysisTree->Branch("ElPt",                   &fTElpt,                "ElPt[NEls]/F");
	fAnalysisTree->Branch("ElEta",                  &fTEleta,               "ElEta[NEls]/F");
	fAnalysisTree->Branch("ElPhi",                  &fTElphi,               "ElPhi[NEls]/F");
	fAnalysisTree->Branch("ElD0",                   &fTEld0,                "ElD0[NEls]/F");
	fAnalysisTree->Branch("ElD0Err",                &fTElD0Err,             "ElD0Err[NEls]/F");
	fAnalysisTree->Branch("ElDz",                   &fTEldz,                "ElDz[NEls]/F");
	fAnalysisTree->Branch("ElDzErr",                &fTElDzErr,             "ElDzErr[NEls]/F");
	fAnalysisTree->Branch("ElDetIso",               &fTElDetIso,            "ElDetIso[NEls]/F");
	fAnalysisTree->Branch("ElPFIso",                &fTElPFIso,             "ElPFIso[NEls]/F");
	fAnalysisTree->Branch("ElPFChIso",              &fTElPFchiso,           "ElPFChIso[NEls]/F");
	fAnalysisTree->Branch("ElPFNeIso",              &fTElPFneiso,           "ElPFNeIso[NEls]/F");
	fAnalysisTree->Branch("ElRadIso",               &fTElRadIso,            "ElRadIso[NEls]/F");
	fAnalysisTree->Branch("ElMVAIDnoTrig",          &fTElMVAIDnoTrig,       "ElMVAIDnoTrig[NEls]/F");
	fAnalysisTree->Branch("ElMVAIDTrig",            &fTElMVAIDTrig,         "ElMVAIDTrig[NEls]/F");
	fAnalysisTree->Branch("ElEcalRecHitSumEt",      &fTElEcalRecHitSumEt,   "ElEcalRecHitSumEt[NEls]/F");
	fAnalysisTree->Branch("ElHcalTowerSumEt",       &fTElHcalTowerSumEt,    "ElHcalTowerSumEt[NEls]/F");
	fAnalysisTree->Branch("ElTkSumPt",              &fTElTkSumPt,           "ElTkSumPt[NEls]/F");
	fAnalysisTree->Branch("ElDPhi",                 &fTElDPhi,              "ElDPhi[NEls]/F");
	fAnalysisTree->Branch("ElDEta",                 &fTElDEta,              "ElDEta[NEls]/F");
	fAnalysisTree->Branch("ElSigmaIetaIeta",        &fTElSigmaIetaIeta,     "ElSigmaIetaIeta[NEls]/F");
	fAnalysisTree->Branch("ElHoverE",               &fTElHoverE,            "ElHoverE[NEls]/F");
	fAnalysisTree->Branch("ElEPthing",              &fTElEPthing,           "ElEPthing[NEls]/F");
	fAnalysisTree->Branch("ElIsGoodElId_LooseWP",      &fTElIsGoodElId_LooseWP,   "ElIsGoodElId_LooseWP[NEls]/I");
	fAnalysisTree->Branch("ElIsGoodElId_MediumWP",     &fTElIsGoodElId_MediumWP,  "ElIsGoodElId_MediumWP[NEls]/I");
	fAnalysisTree->Branch("ElIsGoodTriggerEl",     &fTElIsGoodTriggerEl,  "ElIsGoodTrigger[NEls]/I");
	fAnalysisTree->Branch("ElGenID",                &fTElGenID,             "ElGenID[NEls]/I");
	fAnalysisTree->Branch("ElGenMID",               &fTElGenMID,            "ElGenMID[NEls]/I");
	fAnalysisTree->Branch("ElGenGMID",              &fTElGenGMID,           "ElGenGMID[NEls]/I");
	fAnalysisTree->Branch("ElGenType",              &fTElGenType,           "ElGenType[NEls]/I");
	fAnalysisTree->Branch("ElGenMType",             &fTElGenMType,          "ElGenMType[NEls]/I");
	fAnalysisTree->Branch("ElGenGMType",            &fTElGenGMType,         "ElGenGMType[NEls]/I");
	fAnalysisTree->Branch("ElMT",                   &fTElMT,                "ElMT[NEls]/F");

	// single-tau properties
	fAnalysisTree->Branch("NTaus",                   &fTnqtaus,               "NTaus/I");
	fAnalysisTree->Branch("TauCharge",               &fTTaucharge,            "TauCh[NTaus]/I");
	fAnalysisTree->Branch("TauPt",                   &fTTaupt,                "TauPt[NTaus]/F");
	fAnalysisTree->Branch("TauEta",                  &fTTaueta,               "TauEta[NTaus]/F");
	fAnalysisTree->Branch("TauPhi",                  &fTTauphi,               "TauPhi[NTaus]/F");
	fAnalysisTree->Branch("TauMVAElRej",             &fTTauMVAElRej,          "TauMVAElRej[NTaus]/F");
	fAnalysisTree->Branch("TauTightMuRej",           &fTTauTightMuRej,        "TauTightMuRej[NTaus]/F");
	fAnalysisTree->Branch("TauLCombIsoDB",           &fTTauLCombIsoDB,        "TauLCombIsoDB[NTaus]/F");

	// jet-MET properties
	fAnalysisTree->Branch("pfMET",         &fTpfMET,       "pfMET/F");
	fAnalysisTree->Branch("pfMETPhi",      &fTpfMETphi,    "pfMETPhi/F");
	fAnalysisTree->Branch("pfMETType1",       &fTpfMETType1,     "pfMETType1/F");
	fAnalysisTree->Branch("pfMETType1Phi",    &fTpfMETType1phi,  "pfMETType1Phi/F");
	fAnalysisTree->Branch("NJets",         &fTnqjets,        "NJets/I");
	fAnalysisTree->Branch("JetPt",         &fTJetpt,         "JetPt[NJets]/F");
	fAnalysisTree->Branch("JetEta",        &fTJeteta,        "JetEta[NJets]/F");
	fAnalysisTree->Branch("JetPhi",        &fTJetphi,        "JetPhi[NJets]/F");
	fAnalysisTree->Branch("JetEnergy",     &fTJetenergy,     "JetEnergy[NJets]/F");
	fAnalysisTree->Branch("JetCSVBTag",    &fTJetbtag1,      "JetCSVBTag[NJets]/F");
	fAnalysisTree->Branch("JetProbBTag",   &fTJetbtag2,      "JetProbBTag[NJets]/F");
	fAnalysisTree->Branch("JetArea",       &fTJetArea,       "JetArea[NJets]/F");
	fAnalysisTree->Branch("JetJEC",        &fTJetJEC,        "JetJEC[NJets]/F");
	fAnalysisTree->Branch("JetPartonID",   &fTJetPartonID,   "JetPartonID[NJets]/I");
	fAnalysisTree->Branch("JetPartonFlav", &fTJetPartonFlav, "JetPartonFlav[NJets]/I");
	fAnalysisTree->Branch("JetGenPt",      &fTJetGenpt ,     "JetGenPt[NJets]/F");
	fAnalysisTree->Branch("JetGenEta",     &fTJetGeneta,     "JetGenEta[NJets]/F");
	fAnalysisTree->Branch("JetGenPhi",     &fTJetGenphi,     "JetGenPhi[NJets]/F");
}

void SSDLAnalysis::BookEffTree(){
	fOutputFile->cd();
	fLepEffTree = new TTree("LeptonEfficiency", "LeptonEfficiencyTree");

    // run/sample properties
	fLepEffTree->Branch("Run",              &fLETrun      ,       "Run/I");
	fLepEffTree->Branch("Event",            &fLETevent    ,       "Event/I");
	fLepEffTree->Branch("LumiSec",          &fLETlumi     ,       "LumiSec/I");
	fLepEffTree->Branch("Rho",              &fLETrho      ,       "Rho/F");
	fLepEffTree->Branch("NVrtx",            &fLETnvrtx    ,       "NVrtx/I");
	fLepEffTree->Branch("PUWeight",         &fLETpuweight ,       "PUWeight/F");
	fLepEffTree->Branch("Type",             &fLETtype     ,       "Type/I");
	fLepEffTree->Branch("Pt",               &fLETpt       ,       "Pt/F");
	fLepEffTree->Branch("Eta",              &fLETeta      ,       "Eta/F");
	fLepEffTree->Branch("Phi",              &fLETphi      ,       "Phi/F");
	fLepEffTree->Branch("Iso",              &fLETiso      ,       "Iso/F");
	fLepEffTree->Branch("Pass1",            &fLETpassed1  ,       "Pass1/I");
	fLepEffTree->Branch("Pass2",            &fLETpassed2  ,       "Pass2/I");
	fLepEffTree->Branch("Pass3",            &fLETpassed3  ,       "Pass3/I");
	fLepEffTree->Branch("Pass4",            &fLETpassed4  ,       "Pass4/I");
}

//____________________________________________________________________________
void SSDLAnalysis::Analyze(){
	fHEvCount->Fill(0.);
	FillAnalysisTree();
	if(fDoFillEffTree && !fIsData) FillEffTree();
}
void SSDLAnalysis::FillAnalysisTree(){
	fCounter.fill(fCutnames[0]);
	bool TChiSlepSnu(false);
	bool isRightHanded(false);
	float x(0);
	if (!fIsData){
		// --------- all the necessary histograms for counting nEvents in an SMS scan are put here. it's a bit messy
		// ---------------------------------------------------------------------------------------------------------
		fMsugraCount->Fill(fTR->M0, fTR->M12);
		if (fTR->process > 0 && fTR->process < 11) fProcessCount[(fTR->process)-1]->Fill(fTR->M0, fTR->M12);
		// define some SMS relevant variables here:
		// isRightHanded: chargino1 decays into taus
		// TChiSlepSnu 	: charged leptons come from chargino1 ...
		for (int i = 0; i < fTR->NGenLeptons; i++) {
			if ( abs(fTR->GenLeptonID[i]) == 15 || abs(fTR->GenLeptonID[i]) == 16 ) {  // 15 & 16 is tau family
				if ( abs(fTR->GenLeptonMID[i]) == 1000024 ) isRightHanded = true; // 1000024 == chargino 1
			}

			if ( abs(fTR->GenLeptonID[i]) == 12 || abs(fTR->GenLeptonID[i]) == 14 || abs(fTR->GenLeptonID[i]) == 16) continue; // 12, 14, 16 are neutrinos
			if ( abs(fTR->GenLeptonMID[i]) == 1000024 ) TChiSlepSnu = true; // 1000024 == chargino 1
		}
		x = fTR->MassChi; // this doesn't make a lot of sense, but it's the way it is.
		// MassGlu is the mass of the produced particle
		// MassLSP is in fact the mass of the LSP

		// now finally filling the histograms
		                                   fModelCountAll           -> Fill(fTR->MassGlu, fTR->MassLSP);
		if (!TChiSlepSnu && isRightHanded) fRightHandedSlepCountAll -> Fill(fTR->MassGlu, fTR->MassLSP);
		if (isRightHanded)                 fRightHandedCountAll     -> Fill(fTR->MassGlu, fTR->MassLSP);
		TChiSlepSnu ?                      fTChiSlepSnuCountAll     -> Fill(fTR->MassGlu, fTR->MassLSP) : fTChiSlepSlepCountAll->Fill(fTR->MassGlu, fTR->MassLSP);

		for (int i=0; i<nx; ++i){ 
			if (!TChiSlepSnu  && isRightHanded && x == x_values[i]) fRightHandedSlepCount[i] -> Fill( fTR->MassGlu, fTR->MassLSP);
			if (isRightHanded &&                  x == x_values[i]) fRightHandedCount[i]     -> Fill( fTR->MassGlu, fTR->MassLSP);
			if (TChiSlepSnu   &&                  x == x_values[i]) fTChiSlepSnuCount[i]     -> Fill( fTR->MassGlu, fTR->MassLSP);
			if (!TChiSlepSnu  &&                  x == x_values[i]) fTChiSlepSlepCount[i]    -> Fill( fTR->MassGlu, fTR->MassLSP);
			if (                                  x == x_values[i]) fModelCount[i]           -> Fill( fTR->MassGlu, fTR->MassLSP);
		}

	}
	// initial event selection: good event trigger, good primary vertex...
	if( !IsGoodEvent() ) return;
	fCounter.fill(fCutnames[1]);
	ResetTree();
	
	// Trigger selection
	// if(fIsData && FillTriggers(fHLTPaths) == false) return;
	FillTriggers();
	fCounter.fill(fCutnames[2]);

	// Do object selections
	vector<int> selectedMuInd  = MuonSelection(           &UserAnalysisBase::IsMostBasicMu);
	vector<int> looseMuInd     = MuonSelection(           &UserAnalysisBase::IsLooseMu);
	vector<int> selectedElInd  = ElectronSelection(       &UserAnalysisBase::IsLooseEl);
	vector<int> selectedTauInd = TauSelection(            &UserAnalysisBase::IsLooseTau);
	vector<int> selectedJetInd = PFJetSelection(15., 2.5, &UserAnalysisBase::IsGoodBasicPFJet);
	fTnqmus  = std::min( (int)selectedMuInd .size(), fMaxNmus );
	fTnqels  = std::min( (int)selectedElInd .size(), fMaxNeles);
	fTnqtaus = std::min( (int)selectedTauInd.size(), fMaxNtaus);
	fTnqjets = std::min( (int)selectedJetInd.size(), fMaxNjets);
	int nLooseMus = std::min( (int)looseMuInd .size(), fMaxNmus );

	// Require at least one loose lepton
	// if( (fTnqmus + fTnqels) < 1 ) return;
	if( (nLooseMus + fTnqels) < 1 ) return;
	fCounter.fill(fCutnames[3]);

	// Event and run info
	fTRunNumber   = fTR->Run;
	fTEventNumber = fTR->Event;
	fTLumiSection = fTR->LumiSection;

	if(!fIsData) {
		fTm0   = fTR->M0;
		fTm12  = fTR->M12;
		fTprocess = fTR->process;
		fTmGlu = fTR->MassGlu;
		fTmLSP = fTR->MassLSP;
		TChiSlepSnu   ? fTisTChiSlepSnu = 1 : fTisTChiSlepSnu = 0;
		isRightHanded ? fTisRightHanded = 1 : fTisRightHanded = 0;
	}
	else {
		fTm0      = -1;
		fTm12     = -1;
		fTprocess = -1;
		fTmGlu    = -1;
		fTmLSP    = -1;
		fTisTChiSlepSnu  = -1;
		fTisRightHanded  = -1;
	}
	// Dump basic jet and MET properties
	for(int ind = 0; ind < fTnqjets; ind++){
		int jetindex = selectedJetInd[ind];
		fTJetpt      [ind] = fTR->JPt[jetindex];
		fTJeteta     [ind] = fTR->JEta[jetindex];
		fTJetphi     [ind] = fTR->JPhi[jetindex];
		fTJetenergy  [ind] = fTR->JE[jetindex];
		fTJetbtag1   [ind] = fTR->JnewPFCombinedSecondaryVertexBPFJetTags[jetindex];
		fTJetbtag2   [ind] = fTR->JnewPFJetProbabilityBPFJetTags[jetindex];
		fTJetArea    [ind] = fTR->JArea[jetindex];
		fTJetJEC     [ind] = GetJECUncert(fTR->JPt[jetindex], fTR->JEta[jetindex])/fTR->JEcorr[jetindex];
		fTJetPartonID[ind] = JetPartonMatch(jetindex);
		fTJetPartonFlav[ind] = fTR->JPartonFlavour[jetindex];
		int genjetind = GenJetMatch(jetindex);
		if(genjetind > -1){
			fTJetGenpt [ind] = fTR->GenJetPt [genjetind];
			fTJetGeneta[ind] = fTR->GenJetEta[genjetind];
			fTJetGenphi[ind] = fTR->GenJetPhi[genjetind];
		}
		else{
			fTJetGenpt [ind] = -888.88;
			fTJetGeneta[ind] = -888.88;
			fTJetGenphi[ind] = -888.88;
		}
		
	}

	// Get METs
	fTpfMET     = fTR->PFMET;
	fTpfMETphi  = fTR->PFMETphi;
	fTpfMETType1   = fTR->PFType1MET;
	fTpfMETType1phi  = fTR->PFType1METphi;
	std::pair<float, float> newmet = GetOnTheFlyCorrections();
	// cout << "        leading jet pt, eta from event: " << fTR->JPt[0] << " , " << fTR->JEta[0] << endl;
	// cout << "JChargedMuEnergyFrac[0]: " << fTR->JChargedMuEnergyFrac[0] << endl;
	if (fabs(1-newmet.first/fTpfMETType1) > 0.001) cout << "---------------------------------" << endl;
	if (fabs(1-newmet.first/fTpfMETType1) > 0.001) cout << "Type1MET from event with phi: " << fTpfMETType1 << " " << fTpfMETType1phi << endl;
	if (fabs(1-newmet.first/fTpfMETType1) > 0.001) cout << "Type1MET from marc  with phi: " << newmet.first << " " << newmet.second   << endl;

	// PU correction
	fTrho   = fTR->Rho;
	fTnvrtx = fTR->NVrtx;
	if(!fIsData) fTpuweight = GetPUWeight(fTR->PUnumInteractions);
	else fTpuweight = 1.;

	// Dump muon properties
	for(int i = 0; i < fTnqmus; ++i){
		int index = selectedMuInd[i];
		fTmupt    [i] = fTR->MuPt      [index];
		fTmueta   [i] = fTR->MuEta     [index];
		fTmuphi   [i] = fTR->MuPhi     [index];
		fTmucharge[i] = fTR->MuCharge  [index];
		fTmupfiso [i] = MuPFIso(index);
		fTmupfiso04 [i] = MuPFIso04(index);
		fTmudetiso[i] = fTR->MuRelIso03[index];
		fTmupfchiso[i] = fTR->MuPfIsoR03ChHad[index] / fTR->MuPt[index];
		fTmupfneiso[i] = (fTR->MuPfIsoR03NeHad[index] + fTR->MuPfIsoR03Photon[index] - 0.5*fTR->MuPfIsoR03SumPUPt[index] ) / fTR->MuPt[index];
		fTmupfneisounc[i] = (fTR->MuPfIsoR03NeHad[index] + fTR->MuPfIsoR03Photon[index]) / fTR->MuPt[index];
		fTmuradiso[i] = MuRadIso(index);
		fTmud0    [i] = fTR->MuD0PV    [index];
		fTmudz    [i] = fTR->MuDzPV    [index];
		fTmuptE   [i] = fTR->MuPtE     [index];
		fTmuEMVetoEt [i] = fTR->MuIso03EMVetoEt [index];
		fTmuHadVetoEt[i] = fTR->MuIso03HadVetoEt[index];
		fTmuPassesTightID[i] = IsGoodBasicMu(index);
		
		if(fIsData == false){ // mc truth information
		        fTIsSignalMuon[i] = IsSignalMuon(index, fTmuid[i], fTmumoid[i], fTmugmoid[i])? 1:0;
			pdgparticle mu, mo, gmo;
			GetPDGParticle(mu,  abs(fTmuid     [i]));
			GetPDGParticle(mo,  abs(fTmumoid   [i]));
			GetPDGParticle(gmo, abs(fTmugmoid  [i]));
			fTmutype   [i] = mu.get_type();
			fTmumotype [i] = mo.get_type();
			fTmugmotype[i] = gmo.get_type();
		} else{
			fTmuid     [i] = -888;
			fTmumoid   [i] = -888;
			fTmugmoid  [i] = -888;
			fTmutype   [i] = -888;
			fTmumotype [i] = -888;
			fTmugmotype[i] = -888;
			fTIsSignalMuon[i] = -999;
		}
		
		// Calculate mT:
		TLorentzVector pmu;
		pmu.SetXYZM(fTR->MuPx[index], fTR->MuPy[index], fTR->MuPz[index], 0.105);	
		double ETlept = sqrt(pmu.M2() + pmu.Perp2());
		double METpx  = fTR->PFType1METpx;
		double METpy  = fTR->PFType1METpy;
		fTmuMT[i]     = sqrt( 2*(fTR->PFType1MET*ETlept - pmu.Px()*METpx - pmu.Py()*METpy ));
	}

	// Dump electron properties
	for(int ind = 0; ind < fTnqels; ind++){
		int elindex = selectedElInd[ind];
		fTElcharge          [ind] = fTR->ElCharge                [elindex];
		fTElChargeIsCons    [ind] = fTR->ElCInfoIsGsfCtfScPixCons[elindex];
		fTElpt              [ind] = fTR->ElPt                    [elindex];
		fTEleta             [ind] = fTR->ElEta                   [elindex];
		fTElphi             [ind] = fTR->ElPhi                   [elindex];
		fTEld0              [ind] = fTR->ElD0PV                  [elindex];
		fTElD0Err           [ind] = fTR->ElD0E                   [elindex];
		fTEldz              [ind] = fTR->ElDzPV                  [elindex];
		fTElDzErr           [ind] = fTR->ElDzE                   [elindex];
		fTElPFIso           [ind] = ElPFIso(elindex);
		fTElDetIso          [ind] = relElIso(elindex);
		fTElPFchiso         [ind] = fTR->ElEventelPFIsoValueCharged03PFIdStandard[elindex] / fTR->ElPt[elindex];
		double neutral = fTR->ElEventelPFIsoValueNeutral03PFIdStandard[elindex] + fTR->ElEventelPFIsoValueGamma03PFIdStandard[elindex];
		double rhocorr = fTR->RhoForIso * Aeff(fTR->ElSCEta[elindex]);
		fTElPFneiso         [ind] = TMath::Max(0., neutral - rhocorr) / fTR->ElPt[elindex];
		fTElRadIso          [ind] = ElRadIso(elindex);
		fTElMVAIDnoTrig     [ind] = fTR->ElIDMVANoTrig[elindex];
		fTElMVAIDTrig       [ind] = fTR->ElIDMVATrig[elindex];
		
		if(fIsData == false){ // mc truth information		
		        fTIsSignalElectron[ind] = IsSignalElectron(elindex, fTElGenID  [ind], fTElGenMID [ind], fTElGenGMID[ind])? 1:0;
			pdgparticle el, emo, egmo;
			GetPDGParticle(el,   abs(fTElGenID  [ind]));
			GetPDGParticle(emo,  abs(fTElGenMID [ind]));
			GetPDGParticle(egmo, abs(fTElGenGMID[ind]));
			fTElGenType  [ind] = el.get_type();
			fTElGenMType [ind] = emo.get_type();
			fTElGenGMType[ind] = egmo.get_type();
		}
		else{
			fTElGenID    [ind] = -888;
			fTElGenMID   [ind] = -888;
			fTElGenGMID  [ind] = -888;
			fTElGenType  [ind] = -888;
			fTElGenMType [ind] = -888;
			fTElGenGMType[ind] = -888;   
			fTIsSignalElectron[ind] = -999;
		}
		
		// Calculate mT:
		TLorentzVector pel;
		pel.SetXYZM(fTR->ElPx[elindex], fTR->ElPy[elindex], fTR->ElPz[elindex], 0.0005);
		double ETlept = sqrt(pel.M2() + pel.Perp2());
		double METpx  = fTR->PFType1METpx;
		double METpy  = fTR->PFType1METpy;
		fTElMT[ind]   = sqrt( 2*(fTR->PFType1MET*ETlept - pel.Px()*METpx - pel.Py()*METpy ));

		// Electron ID
		fTElEcalRecHitSumEt[ind] = fTR->ElDR03EcalRecHitSumEt      [elindex];
		fTElHcalTowerSumEt [ind] = fTR->ElDR03HcalTowerSumEt       [elindex];
		fTElTkSumPt        [ind] = fTR->ElDR03TkSumPt              [elindex];
		fTElDPhi           [ind] = fTR->ElDeltaPhiSuperClusterAtVtx[elindex];
		fTElDEta           [ind] = fTR->ElDeltaEtaSuperClusterAtVtx[elindex];
		fTElSigmaIetaIeta  [ind] = fTR->ElSigmaIetaIeta            [elindex];
		fTElHoverE         [ind] = fTR->ElHcalOverEcal             [elindex];
		fTElEPthing          [ind] = fabs(1/fTR->ElCaloEnergy[elindex] - fTR->ElESuperClusterOverP[elindex]/fTR->ElCaloEnergy[elindex]);
		
		fTElIsGoodElId_LooseWP [ind] = IsGoodElId_LooseWP (elindex);
		fTElIsGoodElId_MediumWP[ind] = IsGoodElId_MediumWP(elindex);
		fTElIsGoodTriggerEl    [ind] = IsGoodTriggerEl    (elindex);
	}

	// Dump tau properties
	for(int ind = 0; ind < fTnqtaus; ind++){
		int tauindex = selectedTauInd[ind];
		fTTaucharge     [ind] = fTR->TauCharge                      [tauindex];
		fTTaupt         [ind] = fTR->TauPt                          [tauindex];
		fTTaueta        [ind] = fTR->TauEta                         [tauindex];
		fTTauphi        [ind] = fTR->TauPhi                         [tauindex];
		fTTauMVAElRej   [ind] = fTR->TauElectronMVARejection        [tauindex];
		fTTauTightMuRej [ind] = fTR->TauTightMuonRejection          [tauindex];
		fTTauLCombIsoDB [ind] = fTR->TauLooseCombinedIsoDBSumPtCorr [tauindex];
	}

	fAnalysisTree->Fill();
}
void SSDLAnalysis::FillEffTree(){
	if(fIsData){
		if(fVerbose > 0) cout << "Trying to fill MC lepton efficiencies on data, Stupid..." << endl;
		return;
	}
	ResetEffTree();

	fLETevent    = fTR->Run;
	fLETrun      = fTR->Event;
	fLETlumi     = fTR->LumiSection;
	fLETrho      = fTR->Rho;
	fLETnvrtx    = fTR->NVrtx;
	fLETpuweight = GetPUWeight(fTR->PUnumInteractions);

	// Muon loop
	for(int i = 0; i < fTR->NMus; ++i){
		if(IsSignalMuon(i) == false) continue; // match to signal mu
		if(fTR->MuPt[i] < 5.) continue; // pt cut
		fLETtype   = 0; // mu
		fLETpt     = fTR->MuPt      [i];
		fLETeta    = fTR->MuEta     [i];
		fLETphi    = fTR->MuPhi     [i];
		// fLETiso    = fTR->MuRelIso03[i];
		fLETiso    = MuPFIso(i);
		
		fLETpassed1 = IsTightMuon(1,i)?1:0;
		fLETpassed2 = IsTightMuon(2,i)?1:0;
		fLETpassed3 = IsTightMuon(3,i)?1:0;
		fLETpassed4 = IsTightMuon(4,i)?1:0;
		fLepEffTree->Fill();
	}

	// Electron loop
	for(int i = 0; i < fTR->NEles; ++i){
		if(IsSignalElectron(i) == false) continue; // match to signal el
		if(fTR->ElPt[i] < 5.) continue; // pt cut
		fLETtype   = 1; // mu
		fLETpt     = fTR->ElPt      [i];
		fLETeta    = fTR->ElEta     [i];
		fLETphi    = fTR->ElPhi     [i];
		// fLETiso    = relElIso(i);
		fLETiso    = ElPFIso(i);
		
		fLETpassed1 = IsTightEle(1,i)?1:0;
		fLETpassed2 = IsTightEle(2,i)?1:0;
		fLETpassed3 = IsTightEle(3,i)?1:0;
		fLETpassed4 = IsTightEle(4,i)?1:0;
		fLepEffTree->Fill();
	}
}

//____________________________________________________________________________
void SSDLAnalysis::ResetTree(){
	// sample/run
	fTRunNumber                  = 0;
	fTEventNumber                = 0;
	fTLumiSection                = 0;

	fTm0      = -999.99;
	fTm12     = -999.99;
	fTprocess = -999;
	fTmGlu    = -999.99;
	fTmLSP    = -999.99;
	fTisTChiSlepSnu   = -999.99;
	fTisRightHanded   = -999.99;

	for(size_t i = 0; i < fHLTPathSets.size(); ++i){
		fHLTResults[i]   -2;
		fHLTPrescales[i] -2;
	}
	
	fTrho      = -999.99;
	fTnvrtx    = -999;
	fTpuweight = -999.99;
	
	// muon properties
	fTnqmus = 0;
	for(int i = 0; i < fMaxNmus; i++){
		fTIsSignalMuon  [i] = -999;
		fTmupt          [i] = -999.99;
		fTmueta         [i] = -999.99;
		fTmuphi         [i] = -999.99;
		fTmucharge      [i] = -999;
		fTmupfiso       [i] = -999.99;
		fTmupfiso04     [i] = -999.99;
		fTmudetiso      [i] = -999.99;
		fTmupfchiso     [i] = -999.99;
		fTmupfneiso     [i] = -999.99;
		fTmupfneisounc  [i] = -999.99;
		fTmuradiso      [i] = -999.99;
		fTmud0          [i] = -999.99;
		fTmudz          [i] = -999.99;
		fTmuptE         [i] = -999.99;
		fTmuEMVetoEt    [i] = -999.99;
		fTmuHadVetoEt   [i] = -999.99;
		fTmuPassesTightID[i] = -999.99;
		fTmuid          [i] = -999;
		fTmumoid        [i] = -999;
		fTmugmoid       [i] = -999;
		fTmutype        [i] = -999;
		fTmumotype      [i] = -999;
		fTmugmotype     [i] = -999;
		fTmuMT          [i] = -999.99;
	}

	// electron properties
	fTnqels = 0;
	for(int i = 0; i < fMaxNeles; i++){
		fTIsSignalElectron  [i] = -999;
		fTElcharge          [i] = -999;
		fTElChargeIsCons    [i] = -999;
		fTElpt              [i] = -999.99;
		fTEleta             [i] = -999.99;
		fTElphi             [i] = -999.99;
		fTEld0              [i] = -999.99;
		fTElD0Err           [i] = -999.99;
		fTEldz              [i] = -999.99;
		fTElDzErr           [i] = -999.99;
		fTElDetIso          [i] = -999.99;
		fTElPFIso           [i] = -999.99;
		fTElPFchiso         [i] = -999.99;
		fTElPFneiso         [i] = -999.99;
		fTElRadIso          [i] = -999.99;
		fTElMVAIDnoTrig     [i] = -999.99;
		fTElMVAIDTrig       [i] = -999.99;
		fTElEcalRecHitSumEt [i] = -999.99;
		fTElHcalTowerSumEt  [i] = -999.99;
		fTElTkSumPt         [i] = -999.99;
		fTElDPhi            [i] = -999.99;
		fTElDEta            [i] = -999.99;
		fTElSigmaIetaIeta   [i] = -999.99;
		fTElHoverE          [i] = -999.99;
		fTElEPthing           [i] = -999.99;
		fTElIsGoodElId_LooseWP [i] = -999;
		fTElIsGoodElId_MediumWP[i] = -999;
		fTElIsGoodTriggerEl    [i] = -999;
		fTElMT              [i] = -999.99;
		fTElGenID           [i] = -999;
		fTElGenMID          [i] = -999;
		fTElGenGMID         [i] = -999;
		fTElGenType         [i] = -999;
		fTElGenMType        [i] = -999;
		fTElGenGMType       [i] = -999;
	}

	// tau properties
	fTnqtaus = 0;
	for(int i = 0; i < fMaxNtaus; i++){
		fTTaucharge[i]     = -999;
		fTTaupt[i]         = -999.99;
		fTTaueta[i]        = -999.99;
		fTTauphi[i]        = -999.99;
		fTTauMVAElRej[i]   = -999.99;
		fTTauTightMuRej[i] = -999.99;
		fTTauLCombIsoDB[i] = -999.99;
	}

	// jet-MET properties
	fTnqjets = 0;
	for(int i = 0; i < fMaxNjets; i++){
		fTJetpt [i]       = -999.99;
		fTJeteta[i]       = -999.99;
		fTJetphi[i]       = -999.99;
		fTJetenergy[i]    = -999.99;
		fTJetbtag1[i]     = -999.99;
		fTJetbtag2[i]     = -999.99;
		// MARC fTJetbtag3[i]     = -999.99;
		// MARC fTJetbtag4[i]     = -999.99;
		fTJetArea[i]      = -999.99;
		fTJetJEC[i]       = -999.99;
		fTJetPartonID[i]  = -999;
		fTJetPartonFlav[i]  = -999;
		fTJetGenpt [i]    = -999.99;
		fTJetGeneta[i]    = -999.99;
		fTJetGenphi[i]    = -999.99;
	}
	fTpfMET      = -999.99;
	fTpfMETphi   = -999.99;
	fTpfMETType1      = -999.99;
	fTpfMETType1phi   = -999.99;
}
void SSDLAnalysis::ResetEffTree(){
	fLETevent      = -999;
	fLETrun        = -999;
	fLETlumi       = -999;
	fLETrho        = -999.99;
	fLETnvrtx      = -999;
	fLETpuweight   = -999.99;
	fLETtype       = -999;
	fLETpt         = -999.99;
	fLETeta        = -999.99;
	fLETphi        = -999.99;
	fLETiso        = -999.99;
	fLETpassed1    = -999;
	fLETpassed2    = -999;
	fLETpassed3    = -999;
	fLETpassed4    = -999;
}

//____________________________________________________________________________
// Some same-sign specific object selections
bool SSDLAnalysis::IsSignalMuon(int index){
  int dummy1, dummy2, dummy3;
  return  IsSignalMuon(index,dummy1,dummy2,dummy3);
}
bool SSDLAnalysis::IsSignalMuon(int index, int &muid, int &mumoid, int &mugmoid){

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
	  muid   = fTR->GenLeptonID  [matched];
	  mumoid = fTR->GenLeptonMID [matched];
	  mugmoid= fTR->GenLeptonGMID[matched];
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
	  muid   = fTR->genInfoId[matchedPart];
	  mumoid = mo.get_pdgid() ;
	  mugmoid= gmo.get_pdgid() ;
	  return false;
	}

	//If we really didn't find a match to a lepton or any other particle, set everything to zero
	muid   = 0;
	mumoid = 0;
	mugmoid= 0;
	return false;
}
bool SSDLAnalysis::IsSignalElectron(int index){
  int dummy1, dummy2, dummy3;
  return  IsSignalElectron(index,dummy1,dummy2,dummy3);
}
bool SSDLAnalysis::IsSignalElectron(int index, int &elid, int &elmoid, int &elgmoid){

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
	  elid   = fTR->GenLeptonID  [matched];
	  elmoid = fTR->GenLeptonMID [matched];
	  elgmoid= fTR->GenLeptonGMID[matched];
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
	  elid   = fTR->genInfoId[matchedPart];
	  elmoid = mo.get_pdgid() ;
	  elgmoid= gmo.get_pdgid() ;
	  return false;
	}

	//If we really didn't find a match to a lepton or any other particle, set everything to zero
	elid   = 0;
	elmoid = 0;
	elgmoid= 0;
	return false;	
}
bool SSDLAnalysis::IsTightMuon(int toggle, int index){
	if(toggle == 1){
		if(IsGoodBasicMu(index) == false) return false;
		// MARC if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.1) return false;
		if(MuPFIso(index) > 0.09) return false;
		return true;
	}
	if(toggle == 2){
		if(IsGoodBasicMu(index) == false) return false;
		// MARC if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.1) return false;
		if(corrMuIso(index) > 0.15) return false;
		return true;
	}
	if(toggle == 3){
		if(IsGoodBasicMu(index) == false) return false;
		// MARC if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.1) return false;
		if(MuPFIso(index) > 0.09) return false;
		return true;
	}
	if(toggle == 4){
		if(IsGoodBasicMu(index) == false) return false;
		// MARC if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.1) return false;
		if(corrMuIso(index) > 0.1) return false;
		return true;
	}
	cout << "Choose your toggle!" << endl;
	return false;
}
bool SSDLAnalysis::IsTightEle(int toggle, int index){
	if(toggle == 1){
		if(IsLooseEl(index) == false) return false;
		// MARC if(fTR->ElDR03EcalRecHitSumEt[index]/fTR->ElPt[index] > 0.2) return false;
		if(fTR->ElCInfoIsGsfCtfScPixCons[index] != 1) return false;
		if(IsGoodElId_MediumWP(index) == false) return false;
		if(ElPFIso(index) > 0.15) return false;
		return true;		
	}
	if(toggle == 2){
		if(IsLooseEl(index) == false) return false;
		// MARC if(fTR->ElDR03EcalRecHitSumEt[index]/fTR->ElPt[index] > 0.2) return false;
		if(fTR->ElCInfoIsGsfCtfScPixCons[index] != 1) return false;
		if(IsGoodElId_MediumWP(index) == false) return false;
		if(corrElIso(index) > 0.15) return false;
		return true;
	}
	if(toggle == 3){
		if(IsLooseEl(index) == false) return false;
		// MARC if(fTR->ElDR03EcalRecHitSumEt[index]/fTR->ElPt[index] > 0.2) return false;
		if(fTR->ElCInfoIsGsfCtfScPixCons[index] != 1) return false;
		if(IsGoodElId_MediumWP(index) == false) return false;
		if(ElPFIso(index) > 0.1) return false;
		return true;		
	}
	if(toggle == 4){
		if(IsLooseEl(index) == false) return false;
		// MARC if(fTR->ElDR03EcalRecHitSumEt[index]/fTR->ElPt[index] > 0.2) return false;
		if(fTR->ElCInfoIsGsfCtfScPixCons[index] != 1) return false;
		if(IsGoodElId_MediumWP(index) == false) return false;
		if(corrElIso(index) > 0.1) return false;
		return true;
	}
	cout << "Choose your toggle!" << endl;
	return false;
}

//____________________________________________________________________________
int SSDLAnalysis::JetPartonMatch(int index){
	if(fIsData) return -1;
	// Returns PDG id of matched parton, any of (1,2,3,4,5,21)
	// Unmatched returns 0
	float jpt  = fTR->JPt[index];
	float jeta = fTR->JEta[index];
	float jphi = fTR->JPhi[index];
	int match = 0;
	float mindr = 100;

	////////////////////////////
	////////////////////////////
	// return -2;
	////////////////////////////
	////////////////////////////

	if(fTR->nGenParticles > 1000) return -2;
	for(size_t i = 0; i < fTR->nGenParticles; ++i){
		// Only status 3 particles
		if(fTR->genInfoStatus[i] != 3) continue;
		
		// Restrict to non-top quarks and gluons
		if(abs(fTR->genInfoId[i]) > 5 && fTR->genInfoId[i] != 21) continue;

		// Cutoff for low pt genparticles
		if(fTR->genInfoPt[i]/jpt < 0.1) continue;

		// Minimize DeltaR
		float DR = Util::GetDeltaR(jeta, fTR->genInfoEta[i], jphi, fTR->genInfoPhi[i]);
		if(DR > 0.5) continue;
		if(DR > mindr) continue;
		mindr = DR;
		match = abs(fTR->genInfoId[i]);
	}
	return match;
}
int SSDLAnalysis::GenJetMatch(int index){
	if(fIsData) return -1;
	// Returns index of matched genjet
	// Unmatched returns -1
	float jpt  = fTR->JPt[index];
	float jeta = fTR->JEta[index];
	float jphi = fTR->JPhi[index];

	int match = -1;
	float mindr = 100;

	for(size_t i = 0; i < fTR->NGenJets; ++i){
		// Gen jets are stored with pt > 10
		// Minimize DeltaR
		float DR = Util::GetDeltaR(jeta, fTR->GenJetEta[i], jphi, fTR->GenJetPhi[i]);
		if(DR > 0.5) continue;
		if(DR > mindr) continue;
		mindr = DR;
		match = i;
	}
	return match;
}

//____________________________________________________________________________
// Same-sign specific isolation
double SSDLAnalysis::corrMuIso(int index){
	double newiso = fTR->MuRelIso03[index] - TMath::Log(fTR->MuPt[index]) * (fTR->NVrtx - 1) / (30. * fTR->MuPt[index]);
	if(newiso < 0.) newiso = 0.; // cut off at 0
	return newiso;
}
double SSDLAnalysis::corrElIso(int index){
	double newiso = relElIso(index) - TMath::Log(fTR->ElPt[index]) * (fTR->NVrtx - 1) / (30. * fTR->ElPt[index]);
	if(newiso < 0.) newiso = 0.; // cut off at 0
	return newiso;
}
