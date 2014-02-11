#include "include/FakeAnalyzer.hh"
#include "base/TreeReader.hh"

using namespace std;

FakeAnalyzer::FakeAnalyzer(std::vector<std::string>& fileList, bool isdata, string globaltag) : TreeAnalyzerBase(fileList) {
	fFakeAnalysis = new FakeAnalysis(fTR, isdata, globaltag);
}

FakeAnalyzer::~FakeAnalyzer(){
	delete fFakeAnalysis;
	// MARC if(!fTR->fChain) cout << "FakeAnalyzer ==> No chain!" << endl;
}

// Method for looping over the tree
void FakeAnalyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << nentries << endl;

	if( fMaxEvents > -1 ){
		cout << " will run on " << fMaxEvents << " events..." << endl;
		nentries = fMaxEvents;
	}

	// MARC if( prescale>1 ) cout << " processing only every " << prescale << " events" << endl;
	// MARC for( Long64_t jentry = 0; jentry < nentries; jentry++ ){

    Long64_t jentry=0;
	for ( fTR->ToBegin(); 
			!(fTR->AtEnd()) && (jentry<fMaxEvents || fMaxEvents<0); 
			++(*fTR) ) 
		{
			PrintProgress(jentry++);

			// MARC PrintProgress(jentry);
			// MARC // Prescale processing...
			// MARC if ( prescale>1 && jentry%prescale ) continue;

			fTR->GetEntry(jentry);

			// Upper Pt Hat cut
			if( (fPtHatCut > -1.0) && (fTR->PtHat > fPtHatCut) ) continue;

			// Run processing
			if( fCurRun != fTR->Run ) { // new run
				fCurRun = fTR->Run;
				skipRun = false;
				if ( CheckRun() == false ) skipRun = true;
				else fFakeAnalysis->BeginRun(fCurRun);
				fFakeAnalysis->BeginRun(fCurRun);
			}
			
			// Check if new lumi is in JSON file
			if( fCurLumi != fTR->LumiSection ) { // new lumisection
				fCurLumi = fTR->LumiSection;
				skipLumi = false;
				if ( CheckRunLumi() == false ) skipLumi = true;
			}
			// disable lumi checking in FakeAnalysis if(skipRun || skipLumi) continue;
			fFakeAnalysis->Analyze();
	}
	cout << endl;
}

// Method called before starting the event loop
// there are two function for with/without PU reweighting
void FakeAnalyzer::BeginJob(){
	// fFakeAnalysis->SetOutputDir(fOutputDir);
	fFakeAnalysis->SetOutputFile(fOutputFile);
	fFakeAnalysis->SetVerbose(fVerbose);
	fFakeAnalysis->SetData(fIsData);
	fFakeAnalysis->Begin();
}

void FakeAnalyzer::BeginJob(std::string dataPuFile, std::string mcPuFile){
	// fFakeAnalysis->SetOutputDir(fOutputDir);
	fFakeAnalysis->FakeAnalysis::SetPileUpSrc( dataPuFile , mcPuFile); //HERE
	fFakeAnalysis->SetOutputFile(fOutputFile);
	fFakeAnalysis->SetVerbose(fVerbose);
	fFakeAnalysis->SetData(fIsData);
	fFakeAnalysis->Begin();
}

// Method called after finishing the event loop
void FakeAnalyzer::EndJob(){
	fFakeAnalysis->End();
}
