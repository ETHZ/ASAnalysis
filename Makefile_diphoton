#
# Macro's Makefile
#
# First run 'make deps', then: 'make [-j n] [targets]'
#
# In order to compile JES feature in UserAnalysisBase, run: 'make [-j n] "DOJES=1" [targets]'
#
ROOTCFLAGS     = $(shell root-config --cflags)
ROOTLIBS       = $(shell root-config --libs)
ROOTGLIBS      = $(shell root-config --glibs)
ROOFIT_INCLUDE = $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR  = $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

#INCLUDES       = -I./include -I${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/interface -I$(CMSSW_RELEASE_BASE)/src/
INCLUDES       = -I./include -I${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/interface -I$(CMSSW_RELEASE_BASE)/src/ -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(LHAPATH)/../../../full/include/

CXX            = g++
CXXFLAGS       = -g -fPIC -fno-var-tracking -Wno-deprecated -D_GNU_SOURCE -O2 -std=c++0x $(INCLUDES) 
#CXXFLAGS       = -g -Wno-deprecated -D_GNU_SOURCE -O2 -std=c++0x -ftree-vectorize $(INCLUDES) 
#CXXFLAGS       = -O2 -fno-var-tracking -pedantic -ansi -pthread -pipe -Wno-vla -Werror=overflow -Wstrict-overflow -std=c++0x -msse3 -ftree-vectorize -Wno-strict-overflow -Werror=array-bounds -Werror=format-contains-nul -Werror=type-limits -fvisibility-inlines-hidden -felide-constructors -fmessage-length=0 -ftemplate-depth-300 -Wall -Wno-non-template-friend -Wno-long-long -Wreturn-type -Wunused -Wparentheses -Wno-deprecated -Werror=return-type -Werror=missing-braces -Werror=unused-value -Werror=address -Werror=format -Werror=sign-compare -Werror=write-strings -fdiagnostics-show-option -g -D_GNU_SOURCE -fPIC  $(INCLUDES)
LD             = g++
LDFLAGS        = -g -lgfortran -lz
SOFLAGS        = -O --no_exceptions -shared


CXXFLAGS      += $(ROOTCFLAGS)

LIBS           = $(ROOTLIBS) 

# Juggle with patch versions...
CMSSW_BASE_VERSION = $(filter CMSSW%, $(subst _patch, , $(CMSSW_VERSION) ))

#NGLIBS         = $(ROOTGLIBS) -lMinuit -lMinuit2 -lTreePlayer
NGLIBS         = $(ROOTGLIBS) -L$(ROOFIT_LIBDIR)/ -lMinuit -lMinuit2 -lTreePlayer -lRooFitCore -lRooFit -lTMVA
GLIBS          = $(filter-out -lNew, $(NGLIBS)) 
GLIBS         += -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -L/swshare/cms/slc5_amd64_gcc462/cms/cmssw/$(CMSSW_BASE_VERSION)/lib/$(SCRAM_ARCH) -lFWCoreFWLite -lFWCoreUtilities -lDataFormatsCommon -lDataFormatsFWLite -lCondFormatsJetMETObjects


				#src/helper/Lumi3DReWeighting_standalone.cc 

SRCS           = src/base/TreeClassBase.cc src/base/TreeReader.cc src/base/TreeAnalyzerBase.cc src/base/UserAnalysisBase.cc \
                 src/mcbtagSFuncert.cc \
                 src/helper/AnaClass.cc src/helper/Davismt2.cc src/helper/FPRatios.cc src/helper/PUWeight.C \
                 src/helper/FakeRatios.cc src/helper/MetaTreeClassBase.C src/helper/GoodRunList.C src/helper/BTagSF.cc src/helper/OnTheFlyCorrections.cc \
                 src/helper/TTGammaScaleFactor.cc
                 
OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend

# Rules ====================================
all: RunDiPhotonJetsAnalyzer 

RunDiPhotonJetsAnalyzer: src/exe/RunDiPhotonJetsAnalyzer.C $(OBJS) src/DiPhotonJetsAnalyzer.o src/DiPhotonMiniTree.o
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^ ${LHAPATH}/../../../full/lib/libLHAPDF.a
	mv -f $@ /scratch/$(USER)
	mv /scratch/$(USER)/$@ $@

clean:
	find src -name '*.o' -exec $(RM) -v {} ';' 
	$(RM) RunDiPhotonJetsAnalyzer

purge:
	$(RM) $(OBJS)

deps: $(SRCS) src/DiPhotonJetsAnalyzer.cc src/DiPhotonMiniTree.cc
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

