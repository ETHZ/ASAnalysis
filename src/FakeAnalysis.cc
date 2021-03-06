#include "FakeAnalysis.hh"
#include "helper/Monitor.hh"
#include "helper/PUWeight.h"

#include "JetCorrectionUncertainty.h"

// #include "LHAPDF/LHAPDF.h"


using namespace std;

TString FakeAnalysis::gBaseDir = "/shome/mdunser/workspace/CMSSW_5_3_7_patch5/src/ASAnalysis/";

//____________________________________________________________________________
FakeAnalysis::FakeAnalysis(TreeReader *tr, bool data, string globaltag): UserAnalysisBase(tr, data, globaltag){
	fHLTPaths.clear();

	fGlobalTag = globaltag;

	fCutnames[0] = "All events";
	fCutnames[1] = " ... passes primary Vertex cuts";
	fCutnames[2] = " ... passes triggers (data only)";
	fCutnames[3] = " ... has at least one loose lepton (mc only)";


	fCounter.fill(fCutnames[0], 0.);
	fCounter.fill(fCutnames[1], 0.);
	fCounter.fill(fCutnames[2], 0.);
	fCounter.fill(fCutnames[3], 0.);
}
FakeAnalysis::~FakeAnalysis(){
}

//____________________________________________________________________________
void FakeAnalysis::Begin(const char* filename){
	ReadTriggers(gBaseDir + "FRtriggers.dat");
	ReadPDGTable(gBaseDir + "pdgtable.txt");

	cout << "fakeanalysis ----------------  isdata? "  << (fIsData ? "yes":"no") << endl;
	BookTree();
	fHEvCount = new TH1F("EventCount", "Event Counter", 1, 0., 1.); // count number of generated events
}

//____________________________________________________________________________
void FakeAnalysis::End(){
	fOutputFile->cd();
	fHEvCount->Write();
	fAnalysisTree->Write();
	fOutputFile->Close();
	fCounter.print();
}

//____________________________________________________________________________
void FakeAnalysis::ReadTriggers(const char* triggerfile){
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

		if(!ok) cout << "%% FakeAnalysis::ReadTriggers ==> ERROR: Reading failed!" << endl;
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
void FakeAnalysis::AddTriggerBranches(){
	for(unsigned int i = 0; i < fHLTPathSets.size(); i++){
		HLTPathSet ps = fHLTPathSets[i];
		TString prescalename = ps.name + "_PS";
		if(AddBranch(ps.name.Data(),      "I", &fHLTResults[i])   == false ) exit(-1);
		if(AddBranch(prescalename.Data(), "I", &fHLTPrescales[i]) == false ) exit(-1);
	}
}
bool FakeAnalysis::FillTriggers(){
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
const bool FakeAnalysis::AddBranch( const char* name, const char* type, void* address, const char* size ){
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
void FakeAnalysis::BookTree(){
	fOutputFile->cd();
	fAnalysisTree = new TTree("Analysis", "AnalysisTree");

    // run/sample properties
	fAnalysisTree->Branch("Run"  ,            &fTRunNumber  ,       "Run/I");
	fAnalysisTree->Branch("Lumi" ,            &fTLumiSection,       "Lumi/I");
	fAnalysisTree->Branch("Event",            &fTEventNumber,       "Event/I");

	// HLT triggers
	AddTriggerBranches();
	
	// event properties
	fAnalysisTree->Branch("NVrtx"         ,&fTnvrtx         , "NVrtx/I"                );
	fAnalysisTree->Branch("NTrue"         ,&fTntrue         , "NTrue/I"                );
	fAnalysisTree->Branch("Rho"           ,&fTrho           , "Rho/F"                );
	fAnalysisTree->Branch("PUWeight"      ,&fTpuweight      , "PUWeight/F"             );
	fAnalysisTree->Branch("PUWeightUp"    ,&fTpuweightUp    , "PUWeightUp/F"           );
	fAnalysisTree->Branch("PUWeightDn"    ,&fTpuweightDn    , "PUWeightDn/F"           );
	fAnalysisTree->Branch("GenWeight"     ,&fTGenWeight     , "GenWeight/F"            );

	// single-muon properties
	fAnalysisTree->Branch("MuPt"     , "std::vector<float>", &p_fTmupt          );
	fAnalysisTree->Branch("MuEta"    , "std::vector<float>", &p_fTmueta         );
	fAnalysisTree->Branch("MuPhi"    , "std::vector<float>", &p_fTmuphi         );
	fAnalysisTree->Branch("MuCharge" , "std::vector<int>"  , &p_fTmucharge      );
	fAnalysisTree->Branch("MuD0"     , "std::vector<float>", &p_fTmud0          );
	fAnalysisTree->Branch("MuPFIso"  , "std::vector<float>", &p_fTmupfiso       );
	fAnalysisTree->Branch("MuNeIso"  , "std::vector<float>", &p_fTmuneiso       );
	fAnalysisTree->Branch("MuPhIso"  , "std::vector<float>", &p_fTmuphiso       );
	fAnalysisTree->Branch("MuChIso"  , "std::vector<float>", &p_fTmuchiso       );
	fAnalysisTree->Branch("MuSumPU"  , "std::vector<float>", &p_fTmusumpu       );

	fAnalysisTree->Branch("MuIsVeto"      , "std::vector<bool>", &p_fTmuisveto  );
	fAnalysisTree->Branch("MuIsLoose"     , "std::vector<bool>", &p_fTmuisloose );
	fAnalysisTree->Branch("MuIsTight"     , "std::vector<bool>", &p_fTmuistight );
	fAnalysisTree->Branch("MuIsPrompt"    , "std::vector<bool>", &p_fTmuisprompt);
	fAnalysisTree->Branch("MuID"          , "std::vector<int>" , &p_fTmuid);
	fAnalysisTree->Branch("MuMID"         , "std::vector<int>" , &p_fTmumid);
	fAnalysisTree->Branch("MuGMID"        , "std::vector<int>" , &p_fTmugmid);
	fAnalysisTree->Branch("MuMatchID"     , "std::vector<int>" , &p_fTmumatchid);
	fAnalysisTree->Branch("MuMatchMID"    , "std::vector<int>" , &p_fTmumatchmid);
	fAnalysisTree->Branch("MuMatchGMID"   , "std::vector<int>" , &p_fTmumatchgmid);

	// // single-electron properties
	fAnalysisTree->Branch("ElPt"     , "std::vector<float>", &p_fTelpt          );
	fAnalysisTree->Branch("ElEta"    , "std::vector<float>", &p_fTeleta         );
	fAnalysisTree->Branch("ElPhi"    , "std::vector<float>", &p_fTelphi         );
	fAnalysisTree->Branch("ElCharge" , "std::vector<int>"  , &p_fTelcharge      );
	fAnalysisTree->Branch("ElD0"     , "std::vector<float>", &p_fTeld0          );
	fAnalysisTree->Branch("ElChCo"   , "std::vector<int>"  , &p_fTelchco        );

	fAnalysisTree->Branch("ElPFIso"  , "std::vector<float>", &p_fTelpfiso       );
	fAnalysisTree->Branch("ElNeIso"  , "std::vector<float>", &p_fTelneiso       );
	fAnalysisTree->Branch("ElPhIso"  , "std::vector<float>", &p_fTelphiso       );
	fAnalysisTree->Branch("ElChIso"  , "std::vector<float>", &p_fTelchiso       );

	fAnalysisTree->Branch("ElIsVeto"      , "std::vector<bool>", &p_fTelisveto  );
	fAnalysisTree->Branch("ElIsLoose"     , "std::vector<bool>", &p_fTelisloose );
	fAnalysisTree->Branch("ElIsTight"     , "std::vector<bool>", &p_fTelistight );
	fAnalysisTree->Branch("ElIsPrompt"    , "std::vector<bool>", &p_fTelisprompt);
	fAnalysisTree->Branch("ElID"          , "std::vector<int>" , &p_fTelid);
	fAnalysisTree->Branch("ElMID"         , "std::vector<int>" , &p_fTelmid);
	fAnalysisTree->Branch("ElGMID"        , "std::vector<int>" , &p_fTelgmid);
	fAnalysisTree->Branch("ElMatchID"     , "std::vector<int>" , &p_fTelmatchid);
	fAnalysisTree->Branch("ElMatchMID"    , "std::vector<int>" , &p_fTelmatchmid);
	fAnalysisTree->Branch("ElMatchGMID"   , "std::vector<int>" , &p_fTelmatchgmid);

	// // jet-MET properties
	fAnalysisTree->Branch("pfMET"         ,&fTpfMET         , "pfMET/F"                );
	fAnalysisTree->Branch("pfMETPhi"      ,&fTpfMETphi      , "pfMETPhi/F"             );
	fAnalysisTree->Branch("pfMET1"        ,&fTpfMET1        , "pfMET1/F"                );
	fAnalysisTree->Branch("pfMET1Phi"     ,&fTpfMET1phi     , "pfMET1Phi/F"             );

	fAnalysisTree->Branch("JetPt"         , "std::vector<float>" , &p_fTJetpt         );
	fAnalysisTree->Branch("JetRawPt"      , "std::vector<float>" , &p_fTJetrawpt      );
	fAnalysisTree->Branch("JetEta"        , "std::vector<float>" , &p_fTJeteta        );
	fAnalysisTree->Branch("JetPhi"        , "std::vector<float>" , &p_fTJetphi        );
	fAnalysisTree->Branch("JetEnergy"     , "std::vector<float>" , &p_fTJetenergy     );
	fAnalysisTree->Branch("JetCSVBTag"    , "std::vector<float>" , &p_fTJetCSVtag     );
	fAnalysisTree->Branch("JetPartonFlav" , "std::vector<int>"   , &p_fTJetPartonFlav );
	fAnalysisTree->Branch("JetBetaStar"   , "std::vector<float>" , &p_fTJetBetaStar   );

}

//____________________________________________________________________________
void FakeAnalysis::Analyze(){
	fHEvCount->Fill(0.);
	FillAnalysisTree();
}
void FakeAnalysis::FillAnalysisTree(){

	fCounter.fill(fCutnames[0]);
	// initial event selection: good event trigger, good primary vertex...
	if( !IsGoodEvent() ) return;
	fCounter.fill(fCutnames[1]);
	ResetTree();
	
	// Trigger selection
	if(fIsData && FillTriggers() == false) return;
	FillTriggers();
	fCounter.fill(fCutnames[2]);

	fCounter.fill(fCutnames[3]);

	// Event and run info
	fTRunNumber   = fTR->Run;
	fTEventNumber = fTR->Event;
	fTLumiSection = fTR->LumiSection;

	// Dump basic jet and MET properties
	for(int ind = 0; ind < fTR->NJets; ++ind){
		if(fabs(fTR->JEta[ind]) > 2.5 || fTR->JPt[ind] <  1.) continue;
		if( !IsGoodBasicPFJet(ind) ) continue; // jet selector
		
		float pt = (fGlobalTag == "" ? fTR->JPt[ind]:getNewJetInfo(ind, "pt") ); // new or old pt
		p_fTJetpt      ->push_back( pt) ;
		p_fTJetrawpt   ->push_back( fTR->JPt[ind]/fTR->JEcorr[ind]) ;
		p_fTJeteta     ->push_back( fTR->JEta[ind] );
		p_fTJetphi     ->push_back( fTR->JPhi[ind] );
		p_fTJetenergy  ->push_back( (fGlobalTag == "" ? fTR->JE[ind]:getNewJetInfo(ind, "e") ) ); // new energy
		p_fTJetCSVtag  ->push_back( fTR->JnewPFCombinedSecondaryVertexBPFJetTags[ind] );

		p_fTJetPartonFlav->push_back( fTR->JPartonFlavour[ind]);
		p_fTJetBetaStar  ->push_back( fTR->JBetaStar[ind]     );
		
	}

	// Get METs
	std::pair<float, float> newmet = GetOnTheFlyCorrections();
	fTpfMET1    = fGlobalTag == "" ? fTR->PFType1MET   : newmet.first ;
	fTpfMET1phi = fGlobalTag == "" ? fTR->PFType1METphi: newmet.second;
	fTpfMET     = fTR->PFMET   ;
	fTpfMETphi  = fTR->PFMETphi;

	// fill a met TLorentzVector
	TLorentzVector met;
	met .SetPtEtaPhiM(fTpfMET , 0., fTpfMETphi , 0.);

	fTrho   = fTR->Rho;

	// PU correction
	fTnvrtx = fTR->NVrtx;
	if(!fIsData) {
		fTntrue = fTR->PUnumTrueInteractions;
		fTpuweight   = GetPUWeight    (fTR->PUnumTrueInteractions); // the factor of 1.38 is derived from 20000 ttW events
		fTpuweightUp = GetPUWeightUp  (fTR->PUnumTrueInteractions);
		fTpuweightDn = GetPUWeightDown(fTR->PUnumTrueInteractions);

		fTGenWeight = fTR->GenWeight;
	}
	else {
		fTntrue      = -1;
		fTpuweight   = 1.;
		fTpuweightUp = 1.;
		fTpuweightDn = 1.;
		fTGenWeight  = 1.;
	}

	int nLooseLeptons(0);

	// Dump muon properties
	for(int ind = 0; ind < fTR->NMus; ++ind){
		if(fTR->MuPt[ind] < 10. || fabs(fTR->MuEta[ind]) > 2.4) continue; //save no muons with pt <5, eta > 2.5 or iso > 1.0
	//cout << Form("%i\t%i\t%i\t%.2f\t%.2f\t%.2f\t%i\t%i",fTR->Run, fTR->LumiSection, fTR->Event, fTR->MuPt[ind], fTR->MuEta[ind], fTR->MuPhi[ind], IsLooseMuon(ind), IsTightMuon(ind)) << endl;
		p_fTmupt    ->push_back( fTR->MuPt      [ind] );
		p_fTmueta   ->push_back( fTR->MuEta     [ind] );
		p_fTmuphi   ->push_back( fTR->MuPhi     [ind] );
		p_fTmucharge->push_back( fTR->MuCharge  [ind] );
		p_fTmud0    ->push_back( fTR->MuD0PV    [ind] );
		p_fTmupfiso ->push_back( MuPFIso(ind)         );

		p_fTmuneiso ->push_back( fTR->MuPfIsoR03NeHad  [ind]  );
		p_fTmuphiso ->push_back( fTR->MuPfIsoR03Photon [ind]  );
		p_fTmuchiso ->push_back( fTR->MuPfIsoR03ChHad  [ind]  );
		p_fTmusumpu ->push_back( fTR->MuPfIsoR03SumPUPt[ind]  );

		p_fTmuisveto  ->push_back( IsVetoMuon(ind)  );
		p_fTmuisloose ->push_back( IsLooseMuon(ind) );
		p_fTmuistight ->push_back( IsTightMuon(ind) );
		
		int id(0), mid(0), gmid(0);
		p_fTmuisprompt ->push_back( IsSignalMuon(ind, id, mid, gmid) );
		p_fTmuid       ->push_back(   id );
		p_fTmumid      ->push_back(  mid );
		p_fTmugmid     ->push_back( gmid );
		
		int id1(0), mid1(0), gmid1(0);
		MatchLepton(ind, 0, id1, mid1, gmid1);
		p_fTmumatchid       ->push_back(   id1 );
		p_fTmumatchmid      ->push_back(  mid1 );
		p_fTmumatchgmid     ->push_back( gmid1 );
		
		if(IsLooseMuon(ind)) nLooseLeptons++;

	}


	// Dump electron properties
	for(int ind = 0; ind < fTR->NEles; ind++){
		if(fTR->ElPt[ind] < 10. || fabs(fTR->ElEta[ind]) > 2.5 || ElPFIso(ind) > 1.0 ) continue; //save no electrons with pt <5, eta > 2.5 or iso > 1.0
		
		p_fTelpt     ->push_back( fTR->ElPt    [ind] );
		p_fTeleta    ->push_back( fTR->ElEta   [ind] );
		p_fTelphi    ->push_back( fTR->ElPhi   [ind] );
		p_fTelcharge ->push_back( fTR->ElCharge[ind] );
		p_fTeld0     ->push_back( fTR->ElD0PV  [ind] );
		p_fTelchco   ->push_back( fTR->ElCInfoIsGsfCtfScPixCons[ind] );
		p_fTelpfiso  ->push_back( ElPFIso(ind)       );

		p_fTelneiso  ->push_back( fTR->ElEventelPFIsoValueNeutral03PFIdStandard[ind]    );
		p_fTelphiso  ->push_back( fTR->ElEventelPFIsoValueGamma03PFIdStandard  [ind]    );
		p_fTelchiso  ->push_back( fTR->ElEventelPFIsoValueCharged03PFIdStandard[ind]    );
		
		p_fTelisveto  ->push_back( IsVetoElectron(ind)  );
		p_fTelisloose ->push_back( IsLooseElectron(ind) );
		p_fTelistight ->push_back( IsTightElectron(ind) );

		int id(-1), mid(-1), gmid(-1);
		p_fTelisprompt ->push_back( IsSignalElectron(ind, id, mid, gmid) );
		p_fTelid       ->push_back(   id );
		p_fTelmid      ->push_back(  mid );
		p_fTelgmid     ->push_back( gmid );

		int id1(0), mid1(0), gmid1(0);
		MatchLepton(ind, 1, id1, mid1, gmid1);
		p_fTelmatchid       ->push_back(   id1 );
		p_fTelmatchmid      ->push_back(  mid1 );
		p_fTelmatchgmid     ->push_back( gmid1 );
		
		
		if(IsLooseElectron(ind)) nLooseLeptons++;
	}

	// // Dump Photon properties
	// for(int ind = 0; ind < fTR->NPhotons; ind++){

	// 	if(!IsGoodPhotonEGMLoose(ind)) continue;
	// 	
	// 	p_fTphpt     ->push_back( fTR->PhoPt    [ind] );
	// 	p_fTpheta    ->push_back( fTR->PhoEta   [ind] );
	// 	p_fTphphi    ->push_back( fTR->PhoPhi   [ind] );
	// 	
	// }

	if(nLooseLeptons < 2) return;

	fAnalysisTree->Fill();
}

//____________________________________________________________________________
void FakeAnalysis::ResetTree(){
	// sample/run
	fTRunNumber   = 0;
	fTEventNumber = 0;
	fTLumiSection = 0;

	for(size_t i = 0; i < fHLTPathSets.size(); ++i){
		fHLTResults[i]   -2;
		fHLTPrescales[i] -2;
	}
	
	fTrho        = -999.  ;
	fTnvrtx      = -999   ;
	fTntrue      = -999   ;
	fTpuweight   = -999.99;
	fTpuweightUp = -999.99;
	fTpuweightDn = -999.99;
	fTGenWeight  = -999.99;
	
	// muon properties
	p_fTmupt    = &fTmupt    ; p_fTmupt    ->reserve(fTR->NMus) ; p_fTmupt    ->clear();
	p_fTmueta   = &fTmueta   ; p_fTmueta   ->reserve(fTR->NMus) ; p_fTmueta   ->clear();
	p_fTmuphi   = &fTmuphi   ; p_fTmuphi   ->reserve(fTR->NMus) ; p_fTmuphi   ->clear();
	p_fTmupfiso = &fTmupfiso ; p_fTmupfiso ->reserve(fTR->NMus) ; p_fTmupfiso ->clear();
	p_fTmucharge= &fTmucharge; p_fTmucharge->reserve(fTR->NMus) ; p_fTmucharge->clear();
	p_fTmud0    = &fTmud0    ; p_fTmud0    ->reserve(fTR->NMus) ; p_fTmud0    ->clear();

	p_fTmuneiso = &fTmuneiso ; p_fTmuneiso ->reserve(fTR->NMus) ; p_fTmuneiso ->clear();
	p_fTmuphiso = &fTmuphiso ; p_fTmuphiso ->reserve(fTR->NMus) ; p_fTmuphiso ->clear();
	p_fTmuchiso = &fTmuchiso ; p_fTmuchiso ->reserve(fTR->NMus) ; p_fTmuchiso ->clear();
	p_fTmusumpu = &fTmusumpu ; p_fTmusumpu ->reserve(fTR->NMus) ; p_fTmusumpu ->clear();

	p_fTmuisveto   = &fTmuisveto  ; p_fTmuisveto   ->reserve(fTR->NMus) ; p_fTmuisveto  ->clear();
	p_fTmuisloose  = &fTmuisloose ; p_fTmuisloose  ->reserve(fTR->NMus) ; p_fTmuisloose ->clear();
	p_fTmuistight  = &fTmuistight ; p_fTmuistight  ->reserve(fTR->NMus) ; p_fTmuistight ->clear();
	p_fTmuisprompt  = &fTmuisprompt ; p_fTmuisprompt  ->reserve(fTR->NMus) ; p_fTmuisprompt ->clear();
	p_fTmuid   = &fTmuid ; p_fTmuid  ->reserve(fTR->NMus) ; p_fTmuid ->clear();
	p_fTmumid   = &fTmumid ; p_fTmumid  ->reserve(fTR->NMus) ; p_fTmumid ->clear();
	p_fTmugmid  = &fTmugmid ; p_fTmugmid  ->reserve(fTR->NMus) ; p_fTmugmid ->clear();
	p_fTmumatchid       = &fTmumatchid       ; p_fTmumatchid  ->reserve(fTR->NMus)       ; p_fTmumatchid ->clear()       ;
	p_fTmumatchmid      = &fTmumatchmid      ; p_fTmumatchmid  ->reserve(fTR->NMus)      ; p_fTmumatchmid ->clear()      ;
	p_fTmumatchgmid     = &fTmumatchgmid     ; p_fTmumatchgmid  ->reserve(fTR->NMus)     ; p_fTmumatchgmid ->clear()     ;

	// // electron properties
	p_fTelpt    = &fTelpt    ; p_fTelpt    ->reserve(fTR->NEles); p_fTelpt    ->clear();
	p_fTeleta   = &fTeleta   ; p_fTeleta   ->reserve(fTR->NEles); p_fTeleta   ->clear();
	p_fTelphi   = &fTelphi   ; p_fTelphi   ->reserve(fTR->NEles); p_fTelphi   ->clear();
	p_fTelcharge= &fTelcharge; p_fTelcharge->reserve(fTR->NEles); p_fTelcharge->clear();
	p_fTeld0    = &fTeld0    ; p_fTeld0    ->reserve(fTR->NEles); p_fTeld0    ->clear();
	p_fTelchco  = &fTelchco  ; p_fTelchco  ->reserve(fTR->NEles); p_fTelchco  ->clear();

	p_fTelpfiso = &fTelpfiso ; p_fTelpfiso ->reserve(fTR->NEles); p_fTelpfiso ->clear();
	p_fTelneiso = &fTelneiso ; p_fTelneiso ->reserve(fTR->NEles); p_fTelneiso ->clear();
	p_fTelphiso = &fTelphiso ; p_fTelphiso ->reserve(fTR->NEles); p_fTelphiso ->clear();
	p_fTelchiso = &fTelchiso ; p_fTelchiso ->reserve(fTR->NEles); p_fTelchiso ->clear();

	// // // photon properties
	// p_fTphpt    = &fTphpt    ; p_fTphpt    ->reserve(fTR->NPhotons); p_fTphpt    ->clear();
	// p_fTpheta   = &fTpheta   ; p_fTpheta   ->reserve(fTR->NPhotons); p_fTpheta   ->clear();
	// p_fTphphi   = &fTphphi   ; p_fTphphi   ->reserve(fTR->NPhotons); p_fTphphi   ->clear();

	p_fTelisveto   = &fTelisveto   ; p_fTelisveto   ->reserve(fTR->NEles)  ; p_fTelisveto  ->clear();
	p_fTelisloose  = &fTelisloose  ; p_fTelisloose  ->reserve(fTR->NEles)  ; p_fTelisloose ->clear();
	p_fTelistight  = &fTelistight  ; p_fTelistight  ->reserve(fTR->NEles)  ; p_fTelistight ->clear();
	p_fTelisprompt = &fTelisprompt ; p_fTelisprompt  ->reserve(fTR->NEles) ; p_fTelisprompt ->clear() ;
	p_fTelid       = &fTelid       ; p_fTelid  ->reserve(fTR->NEles)       ; p_fTelid ->clear()       ;
	p_fTelmid      = &fTelmid      ; p_fTelmid  ->reserve(fTR->NEles)      ; p_fTelmid ->clear()      ;
	p_fTelgmid     = &fTelgmid     ; p_fTelgmid  ->reserve(fTR->NEles)     ; p_fTelgmid ->clear()     ;
	p_fTelmatchid       = &fTelmatchid       ; p_fTelmatchid  ->reserve(fTR->NEles)       ; p_fTelmatchid ->clear()       ;
	p_fTelmatchmid      = &fTelmatchmid      ; p_fTelmatchmid  ->reserve(fTR->NEles)      ; p_fTelmatchmid ->clear()      ;
	p_fTelmatchgmid     = &fTelmatchgmid     ; p_fTelmatchgmid  ->reserve(fTR->NEles)     ; p_fTelmatchgmid ->clear()     ;

	// // jet-MET properties
	fTpfMET         = -999.99;
	fTpfMETphi      = -999.99;
	fTpfMET1        = -999.99;
	fTpfMET1phi     = -999.99;

	p_fTJetpt        = &fTJetpt        ; p_fTJetpt        ->reserve(fTR->NJets); p_fTJetpt        ->clear();
	p_fTJetrawpt     = &fTJetrawpt     ; p_fTJetrawpt     ->reserve(fTR->NJets); p_fTJetrawpt     ->clear();
	p_fTJeteta       = &fTJeteta       ; p_fTJeteta       ->reserve(fTR->NJets); p_fTJeteta       ->clear();
	p_fTJetphi       = &fTJetphi       ; p_fTJetphi       ->reserve(fTR->NJets); p_fTJetphi       ->clear();
	p_fTJetenergy    = &fTJetenergy    ; p_fTJetenergy    ->reserve(fTR->NJets); p_fTJetenergy    ->clear();
	p_fTJetCSVtag    = &fTJetCSVtag    ; p_fTJetCSVtag    ->reserve(fTR->NJets); p_fTJetCSVtag    ->clear();
	p_fTJetArea      = &fTJetArea      ; p_fTJetArea      ->reserve(fTR->NJets); p_fTJetArea      ->clear();
	p_fTJetPartonFlav= &fTJetPartonFlav; p_fTJetPartonFlav->reserve(fTR->NJets); p_fTJetPartonFlav->clear();
	p_fTJetBetaStar  = &fTJetBetaStar  ; p_fTJetBetaStar  ->reserve(fTR->NJets); p_fTJetBetaStar  ->clear();

}

// selection function
bool FakeAnalysis::IsVetoElectron(int ind){

	if( fabs(fTR->ElSCEta[ind]) <  1.479 ){ // Barrel
		if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx [ind]) > 0.007)  return false;
		if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx [ind]) > 0.80 )  return false;
		if(fTR->ElSigmaIetaIeta                  [ind]  > 0.01 )  return false;
		if(fTR->ElHcalOverEcal                   [ind]  > 0.15 )  return false;
	}
	if( fabs(fTR->ElSCEta[ind]) >= 1.479 ){ // Endcap
		if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx [ind]) > 0.01 )  return false;
		if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx [ind]) > 0.70 )  return false;
		if(fTR->ElSigmaIetaIeta                  [ind]  > 0.03 )  return false;
		// if(fTR->ElHcalOverEcal                [ind] > 0.07 )   return false;
	}
	if(fabs(fTR->ElDzPV[ind]) > 0.20) return false;
	if(fabs(fTR->ElD0PV[ind]) > 0.04) return false;
	

	return true;
}

bool FakeAnalysis::IsLooseElectron(int ind){
	float epVariable = fabs(1/fTR->ElCaloEnergy[ind] - fTR->ElESuperClusterOverP[ind]/fTR->ElCaloEnergy[ind]);
	if( fabs(fTR->ElSCEta[ind]) <=  1.479 ){ // Barrel
		if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx [ind]) > 0.004) return false;
		if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx [ind]) > 0.06 ) return false;
		if(fTR->ElSigmaIetaIeta                  [ind]  > 0.01 ) return false;
		if(fTR->ElHcalOverEcal                   [ind]  > 0.12 ) return false;
		if(epVariable > 0.05                                   ) return false;
	}
	if( fabs(fTR->ElSCEta[ind]) > 1.479 && fabs(fTR->ElSCEta[ind]) < 2.5 ){ // Endcap
		if(fabs(fTR->ElDeltaEtaSuperClusterAtVtx [ind]) > 0.007) return false;
		if(fabs(fTR->ElDeltaPhiSuperClusterAtVtx [ind]) > 0.03 ) return false;
		if(fTR->ElSigmaIetaIeta                  [ind]  > 0.03 ) return false;
		if(fTR->ElHcalOverEcal                   [ind]  > 0.10 ) return false;
		if(epVariable > 0.05                                   ) return false;
	}
	
	// add the conversion rejection here
	if (!fTR->ElPassConversionVeto[ind]          ) return false;
	if (fTR->ElNumberOfMissingInnerHits[ind] > 0 ) return false;
	if (fTR->ElCInfoIsGsfCtfScPixCons[ind]  != 1 ) return false;
	
	if(fabs(fTR->ElDzPV[ind]) > 0.1) return false;  
	
	//==== here is the part of the selection that I loose:
	//if(fabs(fTR->ElD0PV[ind]) > 0.02) return false; //d0 cut is completely removed
	// fTR->relPfIso[ind]>0.6 return false;   //FIXME! 
	
	return true;
}
bool FakeAnalysis::IsTightElectron(int ind){
  
	if(!IsLooseElectron(ind)) return false;
	
	if(fabs(fTR->ElD0PV[ind]) > 0.02) return false;
	// fTR->relPfIso[ind]>0.15 return false;   //FIXME: 
	
	return true;
}

bool FakeAnalysis::IsVetoMuon(int ind){
	if((fTR->MuIsGlobalMuon[ind] == 0 && fTR->MuIsTrackerMuon[ind] == 0 ))  return false;
	if(fTR->MuIsPFMuon[ind] == 0) return false;
	if(fabs(fTR->MuEta[ind]) > 2.4) return false;
	return true;
}

bool FakeAnalysis::IsLooseMuon(int ind){
	// Basic muon cleaning and ID
	if(fTR->MuIsGlobalMuon[ind] == 0)    return false;
	if(fTR->MuIsPFMuon[ind] == 0)        return false;

	if(fTR->MuNChi2[ind] > 10)           return false;

	if(fTR->MuNGlMuHits.size() > 0){
		if(fTR->MuNGlMuHits [ind] < 1)       return false; // muon.globalTrack()->hitPattern().numberOfValidMuonHits() 
	}
	else{
		if(fTR->MuNMuHits [ind] < 1)         return false; // this is on the outer track: muon.outerTrack()->hitPattern().numberOfValidHits()
	}

    if(fTR->MuNMatchedStations[ind] < 2) return false; // muon.numberOfMatchedStations()

	if(fabs(fTR->MuD0PV[ind])     > 0.2) return false;

	if(fabs(fTR->MuDzPV[ind])     > 0.2) return false;

	if(fTR->MuNPxHits[ind] < 1)          return false; // muon.innerTrack()->hitPattern().numberOfValidPixelHits()

	if(fTR->MuNSiLayers[ind] < 6)        return false; // muon.innerTrack()->hitPattern().trackerLayersWithMeasurement()

	if(MuPFIso(ind) > 1.0) return false;

	// if(fTR->MuNMuHits [ind] < 1)         return false; // this is on the outer track: muon.outerTrack()->hitPattern().numberOfValidHits()
	// not applied if(fTR->MuIso03EMVetoEt[ind]  > 4.00) return false;
	// not applied if(fTR->MuIso03HadVetoEt[ind] > 6.00) return false;

	return true;
}

bool FakeAnalysis::IsTightMuon(int ind){
	if(!IsLooseMuon(ind)) return false;
	if(MuPFIso(ind) > 0.1) return false;
	if(fabs(fTR->MuD0PV[ind]) > 0.01) return false;
	return true;
}

//____________________________________________________________________________
// Some same-sign specific object selections
int FakeAnalysis::findMotherIndex(int ind){
	int tmpm   = fTR->genInfoMo1[ind];
	int tmpmid = fTR->genInfoId[tmpm];
	while(tmpmid == fTR->genInfoId[ind]){
		tmpm   = fTR->genInfoMo1[tmpm];
		tmpmid = fTR->genInfoId[tmpm];
	}
	return tmpm;
}


bool FakeAnalysis::IsGoodPhotonEGMLoose(int i){
	//from https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonID2012 maybe
	if(!IsGoodPhoton(i)) return false;
	float pt = fTR->PhoPt[i];
	float abseta = fabs(fTR->PhoEta[i]);
	float chhadiso  = TMath::Max(fTR->PhoNewIsoPFCharged[i] - fTR->Rho * EffAreaChargedHad(abseta),(float)0.);
	float chneuiso  = TMath::Max(fTR->PhoNewIsoPFNeutral[i] - fTR->Rho * EffAreaNeutralHad(abseta),(float)0.);
	float photoniso = TMath::Max(fTR->PhoNewIsoPFPhoton[i]  - fTR->Rho * EffAreaPhoton(    abseta),(float)0.);
	if(abseta<1.442){//EB
		if(fTR->PhoSigmaIetaIeta[i] > 0.012) return false;
		if(chhadiso  >  2.6                ) return false;
		if(chneuiso  > (3.5 + 0.04  * pt)  ) return false;
		if(photoniso > (1.3 + 0.005 * pt)  ) return false;
	}
	else if(abseta>1.566){//EE
		if(fTR->PhoSigmaIetaIeta[i] > 0.034) return false;
		if(chhadiso  >  2.3                ) return false;
		if(chneuiso  > (2.9 + 0.04  * pt)  ) return false;
		//no photoniso here
	}
	else return false;
	return true;
}

bool FakeAnalysis::IsGoodPhoton(int i){//new
	if( fTR->PhoPt[i] < 20                                        ) return false; // pt cut
	if( fabs(fTR->PhoEta[i])> 2.4                                 ) return false;
	if( fabs(fTR->PhoEta[i])> 1.442 && fabs(fTR->PhoEta[i])<1.566 ) return false; // veto EB-EE gap
	if( fTR->PhoHoverE2012[i] > 0.05                              ) return false;
//	float HoverE2012 = SingleTowerHoverE(i);
//	if( HoverE2012 < -0.5                                         ) return false; // H/E not calculable due missing matched SC
//	if( HoverE2012 > 0.05                                         ) return false; // H/E cut for 2012
	if(!(fTR->PhoPassConversionVeto[i])                           ) return false; // Conversion safe electron veto
	return true;
}

const float FakeAnalysis::EffAreaChargedHad(float abseta){
	abseta=fabs(abseta); // making sure we're looking at |eta|
	if(abseta<1.0)   return 0.012;
	if(abseta<1.479) return 0.010;
	if(abseta<2.0)   return 0.014;
	if(abseta<2.2)   return 0.012;
	if(abseta<2.3)   return 0.016;
	if(abseta<2.4)   return 0.020;
	return 0.012;
}

const float FakeAnalysis::EffAreaNeutralHad(float abseta){
	abseta=fabs(abseta); // making sure we're looking at |eta|
	if(abseta<1.0)   return 0.030;
	if(abseta<1.479) return 0.057;
	if(abseta<2.0)   return 0.039;
	if(abseta<2.2)   return 0.015;
	if(abseta<2.3)   return 0.024;
	if(abseta<2.4)   return 0.039;
	return 0.072;
}

const float FakeAnalysis::EffAreaPhoton(float abseta){
	abseta=fabs(abseta); // making sure we're looking at |eta|
	if(abseta<1.0)   return 0.148;
	if(abseta<1.479) return 0.130;
	if(abseta<2.0)   return 0.112;
	if(abseta<2.2)   return 0.216;
	if(abseta<2.3)   return 0.262;
	if(abseta<2.4)   return 0.260;
	return 0.266;
}

void FakeAnalysis::MatchLepton(int index, int type, int &id, int &moid, int &gmoid){

	float lEta = type == 0 ? fTR->MuEta[index] : fTR->ElEta[index];
	float lPhi = type == 0 ? fTR->MuPhi[index] : fTR->ElPhi[index];
	// look for a generated particle in a given cone
	int matchedPart = -1;
	float mindrPart(100.);
	for(size_t i = 0; i < fTR->nGenParticles; ++i){
		if(fTR->genInfoStatus[i] != 1  && fTR->genInfoStatus[i] != 3 ) continue; // only status 1 or 3. skip the radiations
		float eta = fTR->genInfoEta[i];
		float phi = fTR->genInfoPhi[i];
		float DR = Util::GetDeltaR(eta, lEta, phi, lPhi);
		if(DR > mindrPart) continue; // minimize DR
		mindrPart = DR;
		matchedPart = i;
	}

	bool isGenPart = true;
	bool isSignalGenPart = false;
	if(mindrPart > 0.3) isGenPart = false; // max distance
	if(matchedPart < 0) isGenPart = false; // match unsuccessful

	if (isGenPart){
		// find mother with different ID
		int moIndex = fTR->genInfoMo1[matchedPart];
		while (fTR->genInfoId[matchedPart] == fTR->genInfoId[moIndex]){
			moIndex  = fTR->genInfoMo1[moIndex];
		}
		// find grand mother with different ID than the mother
		int gmoIndex = fTR->genInfoMo1[moIndex];
		while (fTR->genInfoId[moIndex] == fTR->genInfoId[gmoIndex]){
			gmoIndex  = fTR->genInfoMo1[gmoIndex];
		}
		// GetPDGParticle(mo ,  abs(fTR->genInfoId[moIndex] ));
		// GetPDGParticle(gmo,  abs(fTR->genInfoId[gmoIndex]));
		id   = fTR->genInfoId[matchedPart];
		moid = fTR->genInfoId[moIndex] ; //mo.get_pdgid() ;
		gmoid= fTR->genInfoId[gmoIndex]; //gmo.get_pdgid() ;
	}
}


bool FakeAnalysis::IsSignalMuon(int index, int &muid, int &mumoid, int &mugmoid){

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
bool FakeAnalysis::IsSignalElectron(int index, int &elid, int &elmoid, int &elgmoid){

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
