#! /usr/bin/python
import math
import pickle
import os


def ratioWithBinomErrors(numerator, denominator) :
	num = float(numerator)
	den = float(denominator)
	ratio = num / den
	error = math.sqrt(num * (1. - num/den)) / den
	return (ratio, error)
## void SSDLPlotter::ratioWithBinomErrors(float ntight, float nloose, float &ratio, float &error){
## 	ratio = ntight/nloose;
## 	error = TMath::Sqrt( ntight*(1.0-ntight/nloose) ) / nloose;                  // Binomial
## }


def ratioWithPoissErrors(numerator, denominator) :
	num = float(numerator)
	den = float(denominator)
	ratio = num / den
	error = math.sqrt(num*num*(den+num) / (den*den*den))
	return (ratio, error)
## void SSDLPlotter::ratioWithPoissErrors(float ntight, float nloose, float &ratio, float &error){
## 	ratio = ntight/nloose;
## 	error = TMath::Sqrt( ntight*ntight*(nloose+ntight)/(nloose*nloose*nloose) ); // Poissonian	
## }
## void SSDLPlotter::ratioWithAsymmCPErrors(int passed, int total, float &ratio, float &upper, float &lower){
## 	TEfficiency *eff = new TEfficiency("TempEfficiency", "TempEfficiency", 1, 0., 1.);
## 	eff->SetStatisticOption(TEfficiency::kFCP); // Frequentist Clopper Pearson = default
## 	eff->SetConfidenceLevel(0.683); // 1-sigma = default
## 	if( eff->SetTotalEvents(1, total) && eff->SetPassedEvents(1, passed) ){
## 		ratio = eff->GetEfficiency(1);
## 		upper = eff->GetEfficiencyErrorUp(1);
## 		lower = eff->GetEfficiencyErrorLow(1);
## 	}
## 	else{
## 		ratio = 1;
## 		upper = 1;
## 		lower = 0;
## 	};
## 	delete eff;
## }


def save_obj(obj, name):
	with open('obj/' + name + '.pkl', 'w') as file:
		pickle.dump(obj, file, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
	with open('obj/' + name + '.pkl', 'r') as file:
		return pickle.load(file)


def mkdir(path) :
	if not os.path.exists(path) :
		os.mkdir(path)
