#ifndef VALIDATIONPLOTS_HH
#define VALIDATIONPLOTS_HH

#include "helper/AnaClass.hh"

class ValidationPlots{
public:
  ValidationPlots(){};
  ~ValidationPlots(){
    for (unsigned int i=0; i<histo.size(); i++) delete histo[i];
    cleanUp();
  };
  
  std::vector<TH1F*> histo;
  std::vector<TString> name;
  
  int getindex(TString s) { 
    if (histo.size() != name.size()) {
      cout << "SSDLPlotter::ValidationPlots:getindex() ==> ERROR: mismatch in vector size! " << endl;
      return  -1;
    }
    for(unsigned int i=0; i<histo.size(); i++){
      if (s==name[i]) return i;
    }
    
    return -1;
  }
  
  void addhisto(TString n, int nbin, float xmin, float xmax){
    TH1F* h = new TH1F("Val_"+n,"Validation_"+n, nbin, xmin, xmax);
	    histo.push_back(h);
	    name .push_back(n);
  }
  void fillhisto(TString n, float value){	   
    int ind = getindex(n);
    if (ind < 0) {
      cout << "SSDLPlotter::ValidationPlots::fillhisto() ==> ERROR: index " << ind << " not valid!" << endl;
      return;
    }
    histo[ind]->Fill(value);
  }
  void fillhisto(int ind, float value, float weight = 1.){
    if (ind < 0 || ind > histo.size()) {
      cout << "SSDLPlotter::ValidationPlots::fillhisto() ==> ERROR: index " << ind << " not valid!" << endl;
      return;
    }
    histo[ind]->Fill(value);
  }
  void cleanUp(){
    histo.clear();
    name.clear();
  }
};
  
  

#endif
