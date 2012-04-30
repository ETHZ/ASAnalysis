#include "helper/Utilities.hh"
#include "UserAnalysis.hh"

using namespace std;

UserAnalysis::UserAnalysis(TreeReader *tr) : UserAnalysisBase(tr){
    Util::SetStyle();
}

UserAnalysis::~UserAnalysis(){
}

void UserAnalysis::Begin(){
    // Define the output file of histograms
    const char* filename = "histos.root";
    fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	
    // Define the histograms
    fHjetMult = new TH1D("jetMult", "Jet multiplicity", 15, 0, 15);
    fHjetPt   = new TH1D("jetPt", "Pt of jets", 100, 0., 500.);
    fHmuMult  = new TH1D("muMult", "Muon multiplicity", 5, 0, 5);
    fHmuPt    = new TH1D("muPt", "Pt of muons", 100, 0., 300.);
}

void UserAnalysis::Analyze(){
    // Some event selection
    if( fTR->NJets < 1 ) return;		
	
    // Plot some jet quantities
    int nqjets = 0;
    for( int ij = 0; ij < fTR->NJets; ++ij ){
        fHjetPt->Fill(fTR->CAJPt[ij]);

        // Some jet selection
        if( fTR->CAJPt[ij] > 30. ) nqjets++;
    }    
    fHjetMult->Fill(nqjets);

    // Plot some muon quantities
    int nqmus = 0;
    for( size_t im = 0; im < fTR->MuPt.size(); ++im ){
      //if ( !GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v4") ) continue;
      if(fTR->MuIsGlobalMuon[im] == 0) continue;
      fHmuPt->Fill(fTR->MuPt[im]);
      
      // Some muon selection
      if( fTR->MuPt[im] < 10 ) continue;
      if( fTR->MuRelIso03[im] > 1 ) continue;
      if( fTR->MuNTkHits[im] < 11 ) continue;
      if( fTR->MuNChi2[im] > 10 ) continue;
      nqmus++;
    }    
    fHmuMult->Fill(nqmus);

    for ( size_t ih = 0; ih<fTR->HLTObjectPt0.size(); ++ih ) {
      fTR->HLTObjectID0[ih];
    }
}

void UserAnalysis::End(){
    fHistFile->cd();	
    fHjetMult->Write();
    fHjetPt  ->Write();
    fHmuMult ->Write();
    fHmuPt   ->Write();	
    fHistFile->Close();
}
