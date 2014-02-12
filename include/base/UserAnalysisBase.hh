#ifndef UserAnalysisBase_hh
#define UserAnalysisBase_hh

#include <map>
#include <string> 
#include <fstream>

#include <TLatex.h>
#include <TFile.h>
#include <TString.h>

#include "helper/pdgparticle.hh"
#include "helper/Utilities.hh"
#include "helper/LumiReweightingStandAlone.h"

#include "TreeReader.hh"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "helper/OnTheFlyCorrections.hh"

class UserAnalysisBase{
public:
    UserAnalysisBase(TreeReader *tr = 0, bool isData=1, string globaltag="");
    virtual ~UserAnalysisBase();
  
    virtual void Begin() {}
    virtual void BeginRun(Int_t& run);
    virtual void Analyze() {}
    virtual void End() {}
    inline virtual void SetTag(TString tag){fTag = tag;};
    inline virtual void SetVerbose(int verbose){fVerbose = verbose;};
    inline virtual void SetData(bool isdata){fIsData = isdata;};

    inline void SetOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
    inline void SetOutputFile(TString file){ fOutputFile = Util::MakeOutputFile(fOutputDir + file); };

    virtual void ReadPDGTable(const char* filename = "pdgtable.txt");
    virtual int GetPDGParticle(pdgparticle&, int);
    virtual void GetHLTNames();
    virtual int GetHLTBit(string);
    virtual bool GetHLTResult(string);
    virtual int GetHLTPrescale(string);

	OnTheFlyCorrections * fMetCorrector;
    virtual std::pair<float , float > GetOnTheFlyCorrections();

    // PileUp reweighting;
    virtual void  SetPileUpSrc(string, string = "");
    virtual float GetPUWeight(float);
    virtual float GetPUWeightUp(float);
    virtual float GetPUWeightDown(float);

    TreeReader *fTR;
    TString fOutputDir;
    TFile *fOutputFile;
    TString fTag;
    TLatex *fTlat;

    bool fIsData;

    int fVerbose;
    map<int, pdgparticle> fPDGMap; // Mapping of PDG ID names
    map<string, int> fHLTLabelMap; // Mapping of HLT trigger bit names
    vector<string>   fHLTLabels;   // Vector with current HLT names
	
	virtual int   getSusyMass(int, int=1);
	virtual int   getNParticle(int, int=3);
	virtual float getSusySystemPt(int, int=-1);
	virtual float getISRWeight(float, int);

    // Jet Selectors
    virtual bool IsGoodBasicPFJet   (int);
    virtual bool IsGoodPFJetMedium  (int);
    virtual bool IsGoodPFJetTight   (int);

    // Muon Selectors
	virtual float MuPFIso(int);
	virtual float MuPFIso04(int);
	virtual float MuRadIso(int);

    // Electron Selectors
    virtual bool  ElPassesConvRej(int);
    virtual float ElRelDetIso(int);
    virtual float ElPFIso(int);
    virtual float ElRadIso(int);
	virtual float Aeff(float);

    // Event Selectors
    virtual bool IsGoodEvent();

    // Put all JES-related stuff between precompiler flags
	virtual float getNewJetInfo(int ind, string which);
    virtual float GetJetPtNoResidual(int);
	virtual float GetJECUncert(float pt, float eta);

private:

  // Pile UP reweighting
  bool fDoPileUpReweight;
  reweight::LumiReWeighting   *fPUWeight;
  reweight::LumiReWeighting   *fPUWeightUp;
  reweight::LumiReWeighting   *fPUWeightDown;


};

#endif
