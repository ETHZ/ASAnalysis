// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>

#include "TnP.hh"

using namespace std;

//_____________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunTnP [-v verbose] [-i input] [-o output] [-c createHistos]" << endl;
	cout << "  where:" << endl;
	cout << "     verbose         sets the verbose level                  " << endl;
	cout << "     input           is the file containing the histos       " << endl;
	cout << "     output          is the output directory                 " << endl;
	cout << "     createHistos    set to true if creation of histograms is" << endl;
	cout << "                     necessary (takes a while)               " << endl;
	cout << "                     default is false (0)                    " << endl;
	cout << endl;
	exit(status);
}

//_____________________________________________________________________________________
int main(int argc, char* argv[]) {
// Default options
	TString inputfile  = "";
	TString outputdir  = "";
	int verbose = 0;
	bool createHistos = false;

	// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "v:i:o:c:h?")) != -1 ) {
		switch (ch) {
			case 'v': verbose      = atoi(optarg);         break;
			case 'i': inputfile    = TString(optarg);      break;
			case 'o': outputdir    = TString(optarg);      break;
			case 'c': createHistos = (bool) atoi(optarg);  break;
			case '?':
			case 'h': usage(0); break;
			default:
			cerr << "*** Error: unknown option " << optarg << std::endl;
			usage(-1);
		}
	}

	// Check arguments
	if( argc<1 ) {
		usage(-1);
	}
	if(verbose > 0) {
		cout << "------------------------------------" << endl;
		cout << " Verbose level is:  " << verbose << endl;
		cout << " Inputfile is:      " << inputfile << endl;
		cout << " Outputdir is:      " << outputdir << endl;
		if( createHistos) cout << " Create histos: yes   " << endl;
		else              cout << " Create histos: no    " << endl;
	}

	TnP *tA = new TnP(inputfile, createHistos);
	tA->setVerbose(verbose);
	tA->setOutputDir(outputdir);
	tA->doFitting();
	delete tA;
	return 0;
}

