#ifndef TnP_HH
#define Tnp_HH

#include "helper/AnaClass.hh"
#include "helper/Monitor.hh"

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "../../MyRooFitCrap/interface/RooDoubleCB.hh"


using namespace std;

class TnP : public AnaClass{

public:

	TnP(TString, bool);
	virtual ~TnP();
        
	// global variables
	bool fCreateHistos;
	TString fInputFile;
	TString fOutputFileName;

	bool fIsMu;
	bool fIsData;
	int  fnBins;

	// functions
        int Wait();
	virtual void checkFlavor();
	virtual void doFitting();

	virtual void fillHistos();
	virtual void writeHistos();

	virtual void bookHistos();
	virtual void deleteHistos();
	virtual int  readHistos(TString filename);

	TH1F* fMass2LPass;
	TH1F* fMass2LFail;


	virtual void simFitPassFail(TH1F* passHisto , TH1F* failHisto, int flag, int bin);


	float fIsoCut;	
	float fD0Cut;

	virtual bool passesAll(float iso, float d0, float dz, int passid);
	virtual bool passesID(int passid);
	virtual bool passesIP(float d0, float dz);
	virtual bool passesIso(float iso);


	TH1F** fPassIDoverAll  ;
	TH1F** fPassIPoverID   ;
	TH1F** fPassISOoverIDIP;

	TH1F** fFailIDoverAll  ;
	TH1F** fFailIPoverID   ;
	TH1F** fFailISOoverIDIP;

	struct eff{
		int bin;
		float idEff;
		float idEffErr;
		float ipEff;
		float ipEffErr;
		float isoEff;
		float isoEffErr;
	};

	eff* effs;

	// flavor specific things
	// muons
	virtual void printMuTable();
	TString getPtEtaFromIndexMu(int);
	virtual int getPtEtaIndexMu(float pt, float eta);

	// electrons
	virtual void printElTable();
	TString getPtEtaFromIndexEl(int);
	virtual int getPtEtaIndexEl(float pt, float eta);
	virtual bool checkElEta(float);

	virtual void setHistoStyle(TH1F*);
	
private:
	
};

#endif
