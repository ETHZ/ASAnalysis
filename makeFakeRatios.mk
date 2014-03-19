# make shared object of FakeRatios class

python/FakeRatios.so : include/helper/FakeRatios.hh src/helper/FakeRatios.cc LinkdefFakeRatios.h
	rootcint -f FakeRatiosDict.C -c  include/helper/FakeRatios.hh LinkdefFakeRatios.h
	gcc -O -fPIC -I $(ROOTSYS)/include/ -I include/ -c -o FakeRatiosDict.o FakeRatiosDict.C
	gcc -O -fPIC -I $(ROOTSYS)/include/ -I include/ -I /swshare/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/include/ -L /swshare/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms16/lib/ -o FakeRatios.o -c src/helper/FakeRatios.cc
	gcc -shared -lgmp -lgmpxx -o python/FakeRatios.so FakeRatios.o FakeRatiosDict.o
