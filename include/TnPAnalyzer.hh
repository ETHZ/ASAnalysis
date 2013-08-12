#ifndef TnPAnalyzer_hh
#define TnPAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "TnPAnalysis.hh"

class TnPAnalyzer : public TreeAnalyzerBase {
public:
  TnPAnalyzer(std::vector<std::string>& fileList);
  virtual ~TnPAnalyzer();
  void BeginJob(bool, bool);
  void EndJob();
  void Loop();
  void SetOutputFileName(TString a){fTnPAnalysis->SetOutputFileName(a);}



private:
  TnPAnalysis *fTnPAnalysis;


};
#endif
