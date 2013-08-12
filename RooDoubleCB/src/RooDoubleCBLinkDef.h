#pragma GCC diagnostic ignored "-Wwrite-strings" //needed to get rid of pesky "deprecated conversion from string constant to char *" compilation error
#include "../interface/RooDoubleCB.hh"
#include "TVirtualFFT.h"

#ifdef __CINT__

//never even gets here...
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
//#pragma GCC diagnostic ignored "-Wformat"
// #pragma GCC diagnostic warning "-Wwrite-strings"

#pragma link C++ class RooDoubleCB;


#pragma link C++ global gROOT;
#pragma link C++ global gEnv;


#endif
