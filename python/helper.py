#! /usr/bin/python
import math
import pickle
import os
import ROOT
from array import array


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


def getGraphPoissonErrors(histo, nSigma = 1., xErrType = '0') :

	graph = ROOT.TGraphAsymmErrors(0)

	for iBin in range(1, histo.GetXaxis().GetNbins()+1) :

		x = histo.GetBinCenter(iBin)

		if xErrType == '0' :
			xerr = 0.
		elif xErrType == 'binWidth' :
			xerr = histo.GetBinWidth(iBin)/2.
		elif xErrType == 'sqrt12' :
			xerr = histo.GetBinWidth(iBin)/math.sqrt(12.)
		else :
			print '[WARNING] Unkown xErrType %s. Setting to bin width.'
			xerr = histo.GetBinWidth(iBin)

		y = int(histo.GetBinContent(iBin))
		ym = array('d', [0])
		yp = array('d', [0])

		ROOT.RooHistError.instance().getPoissonInterval(y, ym, yp, nSigma)

		yerrplus = yp[0] - y
		yerrminus = y - ym[0]

		thisPoint = graph.GetN()
		graph.SetPoint(thisPoint, x, y )
		graph.SetPointError(thisPoint, xerr, xerr, yerrminus, yerrplus )

	return graph


def save_object(obj, filepath) :
	dir = os.path.dirname(filepath)
	mkdir(dir)
	with open(filepath, 'w') as file :
		pickle.dump(obj, file, pickle.HIGHEST_PROTOCOL)


def load_object(filepath) :
	if not os.path.exists(filepath) :
		print '[ERROR] %s does not exist!' % (filepath)
		return -1
	with open(filepath, 'r') as file :
		return pickle.load(file)


def mkdir(path) :
	if not os.path.exists(path) :
		os.makedirs(path)
