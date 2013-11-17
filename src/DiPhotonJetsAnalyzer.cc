#include "DiPhotonJetsAnalyzer.hh"
#include "DiPhotonMiniTree.hh"


#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"

using namespace std;

DiPhotonJetsAnalyzer::DiPhotonJetsAnalyzer(std::vector<std::string>& fileList, std::string dataType, Float_t aw, Float_t* _kfac, Float_t _minthrpfphotoncandEB, Float_t _minthrpfphotoncandEE, bool _isstep2, TString _input_filename, UInt_t _uuid) : TreeAnalyzerBase(fileList), AddWeight(aw), kfactors(_kfac), minthrpfphotoncandEB(_minthrpfphotoncandEB), minthrpfphotoncandEE(_minthrpfphotoncandEE), isstep2(_isstep2), input_filename(_input_filename), uuid(_uuid) {
  fDiPhotonMiniTree = new DiPhotonMiniTree(fTR,dataType,AddWeight,kfactors,minthrpfphotoncandEB,minthrpfphotoncandEE,isstep2,input_filename,uuid);
}

DiPhotonJetsAnalyzer::~DiPhotonJetsAnalyzer(){
	delete fDiPhotonMiniTree;
}

// Method for looping over the tree
void DiPhotonJetsAnalyzer::Loop(){
    Long64_t nentries = fTR->GetEntries();
    cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	
    // loop over all ntuple entries
    // nentries = 200;
    Long64_t jentry=0;
    for ( fTR->ToBegin(); 
          !(fTR->AtEnd()) && (jentry<fMaxEvents || fMaxEvents<0); 
          ++(*fTR) ) 
        {
            PrintProgress(++jentry);
            if ( fCurRun != fTR->Run ) {
                fCurRun = fTR->Run;
                fDiPhotonMiniTree->BeginRun(fCurRun);
            }
            fDiPhotonMiniTree->Analyze();
        }
}

// Method called before starting the event loop
void DiPhotonJetsAnalyzer::BeginJob(string fdata_PileUp, string fmc_PileUp){

	fDiPhotonMiniTree->SetOutputDir(fOutputDir);
	fDiPhotonMiniTree->SetOutputFile(fOutputFile);
	fDiPhotonMiniTree->fVerbose = fVerbose;
	fDiPhotonMiniTree->SetPileUpSrc(fdata_PileUp, fmc_PileUp);
	fDiPhotonMiniTree->Begin();

}

// Method called after finishing the event loop
void DiPhotonJetsAnalyzer::EndJob(){
	fDiPhotonMiniTree->End();
}

