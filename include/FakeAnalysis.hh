#ifndef FakeAnalysis_hh
#define FakeAnalysis_hh

#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <time.h>

#include <TH2D.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLeaf.h>
#include <TBranch.h>

#include "base/TreeReader.hh"
#include "helper/Utilities.hh"
#include "base/UserAnalysisBase.hh"
#include "helper/Monitor.hh"
#include "helper/PUWeight.h"


using namespace std;

class FakeAnalysis : public UserAnalysisBase{
public:
	// test mar 14 FakeAnalysis(TreeReader *tr = 0);
	FakeAnalysis(TreeReader *tr = 0, bool data =1, string globaltag="");
	virtual ~FakeAnalysis();
	
	void Begin(const char* filename = "SSDLTree.root");
	void Analyze();
	void End();
	
	
	void ReadTriggers(const char* = "HLTPaths_SSDL.dat");
	void AddTriggerBranches();
	bool FillTriggers(); // Returns OR of list of triggers
	void FillAnalysisTree();
	void BookTree();
	void ResetTree();

	int findMotherIndex(int);
	
	const bool AddBranch(const char* name, const char* type, void* address, const char* size = 0);

	struct HLTPathSet{
		TString name;
		vector<string> paths;
	};

	// lepton selection functions
	bool IsVetoElectron (int);
	bool IsLooseElectron(int);
	bool IsTightElectron(int);

	bool IsSignalElectron(int, int&, int&, int&);

	bool IsVetoMuon (int);
	bool IsLooseMuon(int);
	bool IsTightMuon(int);

	bool IsSignalMuon(int, int&, int&, int&);

	// photon function
	bool IsGoodPhoton(int);
	bool IsGoodPhotonEGMLoose(int);

	const float EffAreaChargedHad(float);
	const float EffAreaNeutralHad(float);
	const float EffAreaPhoton(float);

private:
	TH1F *fHEvCount;
	
	
	Monitor fCounter;
	string fCutnames[4];

	static TString gBaseDir; //might be able to replace this at some point
	
	JetCorrectionUncertainty *fJetCorrUnc;
	
	TTree* fAnalysisTree;
	
	std::string fGlobalTag;

	/////////////////////////////////////
	// Tree branches
	int   fTRunNumber;
	int   fTEventNumber;
	int   fTLumiSection;

	// PileUP info
	int   fTnvrtx;
	int   fTntrue;
	float fTpuweight;
	float fTpuweightUp;
	float fTpuweightDn;
	float fTGenWeight;
	
	// Triggers
	vector<HLTPathSet> fHLTPathSets;
	vector<string> fHLTPaths;
	vector<int>    fHLTResults;
	vector<int>    fHLTPrescales;	

	// Muon properties
	vector<float> fTmupt    ; vector<float> * p_fTmupt    ;
	vector<float> fTmueta   ; vector<float> * p_fTmueta   ;
	vector<float> fTmuphi   ; vector<float> * p_fTmuphi   ;
	vector<float> fTmupfiso ; vector<float> * p_fTmupfiso ;
	vector< int > fTmucharge; vector< int > * p_fTmucharge;
	vector<float> fTmud0    ; vector<float> * p_fTmud0    ;

	vector<bool>  fTmuisveto ; vector<bool> * p_fTmuisveto ;
	vector<bool>  fTmuisloose; vector<bool> * p_fTmuisloose;
	vector<bool>  fTmuistight; vector<bool> * p_fTmuistight;

	vector<bool>  fTmuisprompt; vector<bool> * p_fTmuisprompt;
	vector<int>  fTmuid  ; vector<int> * p_fTmuid;
	vector<int>  fTmumid ; vector<int> * p_fTmumid;
	vector<int>  fTmugmid; vector<int> * p_fTmugmid;

	// Electron properties
	vector<float> fTelpt    ; vector<float> * p_fTelpt    ;
	vector<float> fTeleta   ; vector<float> * p_fTeleta   ;
	vector<float> fTelphi   ; vector<float> * p_fTelphi   ;
	vector<float> fTelpfiso ; vector<float> * p_fTelpfiso ;
	vector< int > fTelcharge; vector< int > * p_fTelcharge;
	vector<float> fTeld0    ; vector<float> * p_fTeld0    ;

	vector<bool>  fTelisveto ; vector<bool> * p_fTelisveto ;
	vector<bool>  fTelisloose; vector<bool> * p_fTelisloose;
	vector<bool>  fTelistight; vector<bool> * p_fTelistight;

	vector<bool>  fTelisprompt; vector<bool> * p_fTelisprompt;
	vector<int>  fTelid  ; vector<int> * p_fTelid;
	vector<int>  fTelmid ; vector<int> * p_fTelmid;
	vector<int>  fTelgmid; vector<int> * p_fTelgmid;

	// Electron properties
	vector<float> fTphpt    ; vector<float> * p_fTphpt    ;
	vector<float> fTpheta   ; vector<float> * p_fTpheta   ;
	vector<float> fTphphi   ; vector<float> * p_fTphphi   ;

	// Jet and MET properties
	float fTpfMET;
	float fTpfMETphi;
	float fTpfMET1;
	float fTpfMET1phi;

	vector<float> fTJetpt        ; vector<float> * p_fTJetpt        ;
	vector<float> fTJetrawpt     ; vector<float> * p_fTJetrawpt     ;
	vector<float> fTJeteta       ; vector<float> * p_fTJeteta       ;
	vector<float> fTJetphi       ; vector<float> * p_fTJetphi       ;
	vector<float> fTJetenergy    ; vector<float> * p_fTJetenergy    ;
	vector<float> fTJetCSVtag    ; vector<float> * p_fTJetCSVtag    ;
	vector<float> fTJetArea      ; vector<float> * p_fTJetArea      ;
	vector< int > fTJetPartonFlav; vector< int > * p_fTJetPartonFlav;
	vector<float> fTJetBetaStar  ; vector<float> * p_fTJetBetaStar  ;
	
};
#endif
