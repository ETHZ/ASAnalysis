#include <stdlib.h>
#include <TTree.h>
#include <TString.h>
#include <sstream>
#include <stdio.h>
#include <assert.h>
#include <fstream>

#include "base/TreeAnalyzerBase.hh"

using namespace std;

TreeAnalyzerBase::TreeAnalyzerBase(std::vector<std::string>& fileList) {
    fTR = new TreeReader(fileList);
    fNEntries = fTR->GetEntries();
    fVerbose = 0;
    fMaxEvents = -1;
    fCurRun = -1; // Initialise to dummy value
    fCurLumi = -1;
    skipLumi = false;
    skipRun = false;
}

TreeAnalyzerBase::~TreeAnalyzerBase(){
    delete fTR;
}


// Method that prints the progress at reasonable frequency
void TreeAnalyzerBase::PrintProgress(Long64_t entry){
    int step = 10;
    // Adapt step in powers of 10 (every 10 below 100, every 100 below 1000, etc.)
    Long64_t power = 1;
    for ( size_t i=1; i<10; ++i ) { // up to 10^10...
        power *= 10;
        if ( !(entry/power) ) break;
        step = power;
    }
    if( !(entry%step) ) cout << ">>> Processing event # " << entry << endl;
    else if ( entry==fNEntries-1 )
        cout << ">>> Processing last event # " << entry << endl;
}

// ---------------------------------
// JSON stuff
void TreeAnalyzerBase::ReadJSON(const char* JSONpath) {
    string line;	
    ifstream JSON(JSONpath);
    if(!JSON) cout << "ERROR : JSON FILE DOES NOT EXIST" << endl;
    assert(JSON);
    if (JSON.is_open()){
      string tmp;
      while (JSON.good()){
        getline (JSON,tmp);
	line+=tmp;
      }
    }
    JSON.close();

    line.erase(0,1);
    line.erase(line.size()-1);
	
    size_t index;
    while( (index=line.find("  ")) != string::npos )  line.erase(index,1);

    // split JSON into runs with associated lumis
    string c =", \"";
    vector<string> splitted;
    while(line.find(c) != string::npos){
        int loc = line.find(c);
        string tmp=line.substr(0, (loc));
        splitted.push_back(tmp);	
        line=line.erase(0, tmp.size()+2);
    }	       
    splitted.push_back(line); // for the last entry
	
    // loop over runs (with lumis)
    for(size_t i=0; i<splitted.size(); ++i){
        fCurRunLumi.run =0;
        fCurRunLumi.lumi_min.clear();
        fCurRunLumi.lumi_max.clear();

        // first seperate runs and save them
        if(splitted[i].find(":") != string::npos){
            int loc = splitted[i].find(":");
            string tmp = splitted[i].substr(1,loc-2);
            stringstream ss ( tmp);
            ss >> fCurRunLumi.run;
            splitted[i].erase(0, loc+1);
        }else {cout << "ERROR in ReadJSON: seperate runs" << endl; return;}

        // now we're left with the lumis, remove the surrounding "[]"
        if(splitted[i].find("[[") == string::npos || splitted[i].find("]]") == string::npos){ cout << "ERROR in ReadJSON: lumi parsing" << endl; return;}
        else{
            int loc  = splitted[i].find("[[");
            int loc2 = splitted[i].find("]]");
            splitted[i].erase(loc   , 1);
            splitted[i].erase(loc2  , 1);
        }
        // left with lumis a la [a, b], [c, d], ...
        while(splitted[i].find("[")!= string::npos){
            int loc1 = splitted[i].find("[");
            int loc2 = splitted[i].find("]");
            int loc3 = splitted[i].find(",", loc1);
            string lumi1= splitted[i].substr(loc1+1, loc3-1);
            string lumi2= splitted[i].substr(loc3+1, loc2-1);
            stringstream ss (lumi1); int lumimin;
            ss >> lumimin;	
            fCurRunLumi.lumi_min.push_back(lumimin);
            stringstream ss2 (lumi2); int lumimax;
            ss2 >> lumimax;
            fCurRunLumi.lumi_max.push_back(lumimax);
            splitted[i].erase(loc1, loc2-loc1+3);
        }
	// push_back run and lumis
	fRunLumis.push_back(fCurRunLumi);
    }

    if(fVerbose >1){
        cout << "-----------------------------------------------" << endl;
        cout << "Reading the following runs and lumi sections from JSON file:" << endl; 
        for(size_t i=0; i<fRunLumis.size(); ++i){
            cout << "  Run: " << fRunLumis[i].run << " Lumis: " ;
            for(size_t k=0; k<fRunLumis[i].lumi_min.size(); ++k){
                cout << " [" << fRunLumis[i].lumi_min[k] << ", " << fRunLumis[i].lumi_max[k] << "] ";
            }
            cout << endl;
        }
        cout << "-----------------------------------------------" << endl;
    }
}

const bool TreeAnalyzerBase::CheckRunLumi(void) const {
    bool good(false);
    if(fRunLumis.size()==0) {
        if (fVerbose>1) 
            cout << "CheckRunLumi: fRunLumis is empty. Assuming no JSON." << endl; 
        return true;
    }
    for(size_t r=0; r<fRunLumis.size(); ++r){
        if(fTR->Run != fRunLumis[r].run) continue;
        for(size_t l=0; l<fRunLumis[r].lumi_min.size(); ++l){
            if(fTR->LumiSection < fRunLumis[r].lumi_min[l]) continue;
            if(fTR->LumiSection > fRunLumis[r].lumi_max[l]) continue;
            good = true;
        }
    }
    if(fVerbose > 3){
        if(good) cout <<  "JSON: accepted ";
        else     cout <<  "JSON: rejected ";
        cout <<  fTR->Run << ":" << fTR->LumiSection << ":" << fTR->Event << endl;
    }
    return good;
}

//________________________________________________________________________________________
const bool TreeAnalyzerBase::CheckRun(void) const {
    // Check if current run is in JSON file
    bool runInJson(false);

    // Not initialized: assume no JSON
    if(fRunLumis.size()==0) {
        if (fVerbose>1) 
            cout << "CheckRunLumi: fRunLumis is empty. Assuming no JSON." << endl; 
        return true;
    }
    for(size_t r=0; r<fRunLumis.size(); ++r)
        if(fTR->Run == fRunLumis[r].run)
            runInJson = true;

    if ( !runInJson && fVerbose>2 ) 
        cout << "Run " << fTR->Run << " not in JSON file" << endl;
    return runInJson;

}
