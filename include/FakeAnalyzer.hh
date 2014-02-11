#ifndef FakeAnalyzer_hh
#define FakeAnalyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "FakeAnalysis.hh"

class FakeAnalyzer : public TreeAnalyzerBase {
public:
	FakeAnalyzer(std::vector<std::string>& fileList, bool isdata, string globaltag="");
	virtual ~FakeAnalyzer();
	void BeginJob	();
	void BeginJob	( std::string, std::string );
	void EndJob		();
	void Loop		();

	inline void SetPtHatCut(float cut){fPtHatCut = cut;};
	inline void DoFillEffTree(bool fill){fDoFillEffTree = fill;};

private:
	FakeAnalysis	*fFakeAnalysis;
	float fPtHatCut;
	bool fDoFillEffTree;
};
#endif
