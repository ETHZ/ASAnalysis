#include "DiPhotonJetsAnalyzer.hh"
#include "DiPhotonPurity.hh"


#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

DiPhotonJetsAnalyzer::DiPhotonJetsAnalyzer(TTree *tree, std::string dataType) : TreeAnalyzerBase(tree) {
  fDiPhotonPurity = new DiPhotonPurity(fTR,dataType);
}

DiPhotonJetsAnalyzer::~DiPhotonJetsAnalyzer(){
	delete fDiPhotonPurity;
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
		  fDiPhotonPurity->BeginRun(fCurRun);
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
		  fDiPhotonPurity->Analyze();
		}
	}
}

// Method called before starting the event loop
void DiPhotonJetsAnalyzer::BeginJob(string fdata_PileUp, string fmc_PileUp){
	fDiPhotonPurity->SetOutputDir(fOutputDir);
	fDiPhotonPurity->fVerbose = fVerbose;
	fDiPhotonPurity->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
	fDiPhotonPurity->Begin();

}

// Method called after finishing the event loop
void DiPhotonJetsAnalyzer::EndJob(){
	fDiPhotonPurity->End();
}