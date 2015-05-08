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


def ratio_withError(numerator, numerator_err, denominator, denominator_err) :
	num     = float(numerator      )
	num_err = float(numerator_err  )
	den     = float(denominator    )
	den_err = float(denominator_err)
	ratio     = num / den
	ratio_err = ratio * (num_err/num) / (den_err/den)
	return (ratio, ratio_err)


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


def getGraphPoissonErrors_new(histo, x_errors = False) :
	'''Asymmetric Error Bars for Poisson Event Counts'''

	alpha = 1 - 0.6827
	graph = ROOT.TGraphAsymmErrors(histo)
	for i in range(0, graph.GetN()) :
		N = graph.GetY()[i]
		L = 0
		if N > 0 : L = ROOT.Math.gamma_quantile(alpha/2,N,1.)
		U =  ROOT.Math.gamma_quantile_c(alpha/2,N+1,1)
		graph.SetPointEYlow(i, N-L)
		graph.SetPointEYhigh(i, U-N)
		if not x_errors :
			graph.SetPointEXlow(i, 0.)
			graph.SetPointEXhigh(i, 0.)

	return graph


def get_signifErrDigits(err) :
	err = abs(err)
	err_str = '%e' % err
	coeff = float(err_str.split('e')[0])
	if coeff < 3.55 :
		return 2
	else :
		return 1


def get_precision(num, err) :
	err_precision = get_signifErrDigits(err)
	precision = get_exponent(num) - get_exponent(err) + get_signifErrDigits(err)
	return [precision, err_precision]


def get_exponent(num) :
	num_str = '%e' % num
	exp = int(num_str.split('e')[1])
	return exp


def get_roundedNumber(num, err, rnd = None, float_digits = None) :
	if rnd == None :
		rnd = get_signifErrDigits(err) - get_exponent(err) - 1
	if float_digits == None :
		float_digits = rnd
	if float_digits < 0 : float_digits = 0
	num_str = format(round(num, rnd), '.%df' % float_digits)
	err_str = format(round(err, rnd), '.%df' % float_digits)
	return (num_str, err_str)


def get_roundedNumberErrors(num, errors) :
	'''returns tuple of strings with rounded value and uncertainties'''
	rnds = []
	err_str = []
	for err in errors :
		rnd = get_signifErrDigits(err) - get_exponent(err) - 1
		rnds.append(rnd)
		if rnd < 0 : float_digits = 0
		else       : float_digits = rnd
		err_str.append(format(round(err, rnd), '.%df' % float_digits))
	rnd = max(rnds)
	if rnd < 0 : float_digits = 0
	else       : float_digits = rnd
	num_str = format(round(num, rnd), '.%df' % float_digits)
	return (num_str,) + tuple(err_str)


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
