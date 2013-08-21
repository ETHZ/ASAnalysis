CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src/pandolf/CommonTools -I$(CMSSW_BASE)/src/ASAnalysis/include
#INCLUDES       = -I./include -I${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/interface -I$(CMSSW_RELEASE_BASE)/src/ -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(LHAPATH)/../../../full/include/

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -lFWCoreFWLite -lCondFormatsJetMETObjects -lTreePlayer
#EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -lTreePlayer




scanZBi: scanZBi.cpp DrawBase.o fitTools.o ../src/SSDLPlotter.o ../src/SSDLDumper.o ZBiCalculator.o
	$(CC) -Wall $(INCLUDES) -o scanZBi scanZBi.cpp *.o $(ROOTFLAG) $(EXTRALIBS)

computeTTWZCrossSection: computeTTWZCrossSection.cpp DrawBase.o fitTools.o ../src/SSDLPlotter.o ../src/SSDLDumper.o ZBiCalculator.o
	$(CC) -Wall $(INCLUDES) -o computeTTWZCrossSection computeTTWZCrossSection.cpp *.o $(ROOTFLAG) $(EXTRALIBS)



DrawBase.o: $(CMSSW_BASE)/src/pandolf/CommonTools/DrawBase.C fitTools.o
	$(CC) $(CFLAGS) $(INCLUDES) fitTools.o $(CMSSW_BASE)/src/pandolf/CommonTools/DrawBase.C $(ROOTFLAG) $(EXTRALIBS) 

fitTools.o: $(CMSSW_BASE)/src/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/pandolf/CommonTools/fitTools.C $(ROOTFLAG) $(EXTRALIBS)

StatTools.o: $(CMSSW_BASE)/src/pandolf/CommonTools/StatTools.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/pandolf/CommonTools/StatTools.cc $(ROOTFLAG) $(EXTRALIBS)

ZBiCalculator.o: ZBiCalculator.cc
	$(CC) $(CFLAGS) $(INCLUDES) ZBiCalculator.cc $(ROOTFLAG) $(EXTRALIBS)

#SSDLPlotter.o: ../src/SSDLPlotter.cc
#	$(CC) $(CFLAGS) $(INCLUDES) ../src/SSDLPlotter.cc $(ROOTFLAG) $(EXTRALIBS)
#
#SSDLDumper.o: ../src/SSDLDumper.cc
#	$(CC) $(CFLAGS) $(INCLUDES) ../src/SSDLDumper.cc $(ROOTFLAG) $(EXTRALIBS)

