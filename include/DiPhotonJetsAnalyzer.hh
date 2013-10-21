#ifndef DiPhotonJetsAnalyzer_hh
#define DiPhotonJetsAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "DiPhotonMiniTree.hh"
#include <string>

class DiPhotonJetsAnalyzer : public TreeAnalyzerBase {
public:
  DiPhotonJetsAnalyzer(std::vector<std::string>& fileList, std::string dataType="data", Float_t aw=-999, Float_t* _kfac=NULL, Float_t minthrpfphotoncandEB=0, Float_t minthrpfphotoncandEE=0, bool _isstep2 = false, TString _input_filename = "", UInt_t _uuid = 0);
	virtual ~DiPhotonJetsAnalyzer();
	void BeginJob(string data_PileUp, string mc_PileUp);
	void EndJob();
	void Loop();
  void SetMaxEvents(int a){fMaxEvents=a;}

private:
  DiPhotonMiniTree *fDiPhotonMiniTree;
  int fMaxEvents;
  Float_t AddWeight;
  Float_t* kfactors;
  std::string templateChoice;
  Float_t minthrpfphotoncandEB;
  Float_t minthrpfphotoncandEE;
  bool isstep2;
  TString input_filename;
  UInt_t uuid;
};
#endif
