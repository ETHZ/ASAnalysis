#ifndef TnP_HH
#define Tnp_HH

#include "helper/AnaClass.hh"
#include "helper/Monitor.hh"

#include "TRandom3.h"
#include "TLorentzVector.h"


using namespace std;

class TnP : public AnaClass{

public:

	TnP(TString, bool);
	virtual ~TnP();
        
	// global variables
	bool fCreateHistos;
	TString fInputFile;
	TString fOutputFileName;

	// functions
	virtual void doFitting();

	virtual void fillHistos();
	virtual void writeHistos();

	virtual void bookHistos();
	virtual void deleteHistos();
	virtual int  readHistos(TString filename);

	TH1D* fMass2LPass;
	TH1D* fMass2LFail;


	virtual void simFitPassFail(TH1D* passHisto , TH1D* failHisto);


	// muon stuff
	float fMuIso;	
	float fMuD0;
	virtual bool passesAllMu(float iso, float d0, float dz, int passid);
	virtual bool passesIDMu(int passid);
	virtual bool passesIPMu(float d0, float dz);
	virtual bool passesIsoMu(float iso);

	virtual int getPtEtaIndexMu(float pt, float eta);

	// histograms to be filled
	// 6 bins in pt and 2 bins in eta. the last one is inclusive in pt and eta
	const static int fnBinsMu = 13;
	TH1D* fMuPassIDoverAll   [fnBinsMu];
	TH1D* fMuPassIPoverID    [fnBinsMu];
	TH1D* fMuPassISOoverIDIP [fnBinsMu];

	TH1D* fMuFailIDoverAll   [fnBinsMu];
	TH1D* fMuFailIPoverID    [fnBinsMu];
	TH1D* fMuFailISOoverIDIP [fnBinsMu];

	struct muEff{
		int bin;
		float idEff;
		float idEffErr;
		float ipEff;
		float ipEffErr;
		float isoEff;
		float isoEffErr;
	};

	muEff muEffs[fnBinsMu];

	// electron stuff

	
private:
	
};

#endif
