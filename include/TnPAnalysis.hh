#ifndef UserAnalysis_hh
#define UserAnalysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1D.h>
#include <TTree.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"

class TnPAnalysis : public UserAnalysisBase{
public:
	TnPAnalysis(TreeReader *tr = NULL);
	virtual ~TnPAnalysis();
	

	void Begin();
	void Analyze();
	void End();
	
	
	void BookTree();
	void FillAnalysisTree(int indexLepTag, int indexLepProbe, float massLep);
	
	void SetOutputFileName(TString a){fOutputFileName=a;}

	void SetIsMu(bool);
	bool fIsMu;

	// void SetMyIsData(bool);
	// bool fMIsData;

	bool LepPassesID(int, bool);
	
	void LookForTagAndProbe(); // this function does all the work actually

	bool PassOnlyIDMu(int);
	bool PassOnlyIDEl(int);

	float QuickHT();

	bool IsSignalElectron(int, float&, float&, float&);
	bool IsSignalMuon(int, float&, float&, float&);
	float getGenZMass();

	// class lepton{
	// 	float pt;
	// 	float eta;
	// 	float phi;
	// 	float mass;
	// 	float iso;
	// 	float d0;
	// 	float dz;
	// 	int passID;
	// 	TLorentzVector p4;
	// 	public:
	// 		setp4() { p4.SetPtEtaPhiM(pt, eta, phi, mass); }
	// };

	private:
	  
	// file for histograms:
	TFile *fHistFile;
	TString fOutputFileName;
	
	TH1D *fHCounter;
	TH1D *fHPtHlt1;
	TH1D *fHEtaHlt1;
	TH1D *fHPhiHlt1;
	
	TH1D *fHPtHlt2;
	TH1D *fHEtaHlt2;
	TH1D *fHPhiHlt2;
	
	TH1D *fHDeltaR;
	TH1D *fHPtTag;
	TH1D *fHMass2LAll;   TH1D *fHMass2LPass;   TH1D *fHMass2LFail;
	
	
	TTree* fAnalysisTree;
	/////////////////////////////////////
	// Tree branches
	int   fTRunNumber;
	int   fTEventNumber;
	int   fTLumiSection;
	int   fTisData;
	int   fTisMuEvent;
	
	// tag variables
	float fTagPt;
	float fTagEta;
	float fTagPhi;
	float fTagIsoRel;
	float fTagD0;
	float fTagDz;
	float fTag3DIP;
	float fTagPassID;
	// probe variables
	float fProbePt;
	float fProbeEta;
	float fProbePhi;
	float fProbePtErr;
	float fProbeIsoRel;
	float fProbeD0;
	float fProbeDz;
	float fProbe3DIP;
	float fProbePassID;
	int   fProbeIsGlobal;
	// other variables
	float fMass2L;
	float fMass2LMatchedLeptons;
	float fMass2LFromZ;
	float fDeltaR;
	float fNvtx ;
	float fRho;
	float fPfMet;
	float fHt;
	int fNTrueInt;
};
#endif
