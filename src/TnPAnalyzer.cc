#include "TnPAnalyzer.hh"
#include "TnPAnalysis.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

TnPAnalyzer::TnPAnalyzer(std::vector<std::string>& fileList) 
    : TreeAnalyzerBase(fileList) {
    fTnPAnalysis = new TnPAnalysis(fTR);
}

TnPAnalyzer::~TnPAnalyzer(){
    delete fTnPAnalysis;
}

// Method for looping over the tree
void TnPAnalyzer::Loop(){
    Long64_t nentries = fTR->GetEntries();
    cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	
    // loop over all ntuple entries
    // nentries = 200;
    Long64_t jentry=0;
    for ( fTR->ToBegin(); 
          !(fTR->AtEnd()) && (jentry<fMaxEvents || fMaxEvents<0); 
          ++(*fTR) ) 
        {
            PrintProgress(jentry++);
            if ( fCurRun != fTR->Run ) {
                fCurRun = fTR->Run;
                fTnPAnalysis->BeginRun(fCurRun);
            }
            fTnPAnalysis->Analyze();
        }
}

// Method called before starting the event loop
void TnPAnalyzer::BeginJob(bool isMu, bool isData){
    fTnPAnalysis->SetOutputDir(fOutputDir);
    fTnPAnalysis->SetVerbose(fVerbose);

    fTnPAnalysis->Begin();
	fTnPAnalysis->SetIsMu(isMu);
	fTnPAnalysis->SetData(isData);
	//fTnPAnalysis->SetMyIsData(isData);

}

// Method called after finishing the event loop
void TnPAnalyzer::EndJob(){
    fTnPAnalysis->End();
}
