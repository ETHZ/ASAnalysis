#! /usr/bin/python
import math
import pickle
import os
import ROOT
from array import array
import tables
import numpy as np
import copy


def ratioWithBinomErrors(numerator, denominator) :
	num = float(numerator)
	den = float(denominator)
	ratio = num / den
	error = math.sqrt(num * (1. - num/den)) / den
	return (ratio, error)


def ratioWithPoissErrors(numerator, denominator) :
	num = float(numerator)
	den = float(denominator)
	ratio = num / den
	error = math.sqrt(num*num*(den+num) / (den*den*den))
	return (ratio, error)


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


def product_withError(factor1, factor1_err, factor2, factor2_err) :
	fact1     = float(factor1    )
	fact1_err = float(factor1_err)
	fact2     = float(factor2    )
	fact2_err = float(factor2_err)
	prod      = fact1 * fact2
	prod_err  = prod * math.sqrt((fact1_err/fact1)**2 + (fact2_err/fact2)**2)
	return (prod, prod_err)


def ratio_withError(numerator, numerator_err, denominator, denominator_err) :
	num     = float(numerator      )
	num_err = float(numerator_err  )
	den     = float(denominator    )
	den_err = float(denominator_err)
	ratio     = num / den
	ratio_err = ratio * math.sqrt((num_err/num)**2 + (den_err/den)**2)
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


def get_deltaPhi(phi1, phi2) :
	result = phi1 - phi2
	while result >   math.pi : result -= 2*math.pi
	while result <= -math.pi : result += 2*math.pi
	return abs(result)


def get_deltaR(eta1, eta2, phi1, phi2) :
	deta = eta1 - eta2
	dphi = get_deltaPhi(phi1, phi2)
	return math.sqrt(deta**2 + dphi**2)


def save_histo2table(histos, processes, path, var = 'bincentre', lumi = '', bin_width = True, last_bin = True, asymmErr = False, x_range = None) :
	filepath = copy.deepcopy(path)
	path = '/'.join(filepath.split('/')[:-1])
	filename = filepath.split('/')[-1]
	if filename.endswith('.dat') :
		filename = filename.rstrip('dat').rstrip('.')
	mkdir(path)
	if lumi != '' :
		lumi_str   = '\t%4.1f' % lumi
		lumi_title = '\tlumi'
	else :
		lumi_str   = ''
		lumi_title = ''
	if asymmErr != False :
		histos_tmp = {}
		for process in asymmErr :
			histos_tmp[process] = histos[process].Clone()
			histos_tmp[process].SetBinErrorOption(ROOT.TH1.kPoisson)
	nbins = histos[processes[0]].GetNbinsX()
	with open(filepath, 'w') as file :
		file.write('binlow\t%s\t%s\t%s_err' % (var, '\t'.join(processes), '_err\t'.join(processes)))
		if asymmErr != False :
			file.write('\t%s_uerr\t%s_derr' % ('\t_uerr'.join(asymmErr), '_derr\t'.join(asymmErr)))
		if bin_width : file.write('\tbinwidth')
		file.write('%s\n' % lumi_title)
		for bin in range(1, nbins+1) :
			file.write('%f\t%f' % (histos[processes[0]].GetXaxis().GetBinLowEdge(bin), histos[processes[0]].GetXaxis().GetBinCenter(bin)))
			for process in processes :
				file.write('\t%f' % histos[process].GetBinContent(bin))
			for process in processes :
				file.write('\t%f' % histos[process].GetBinError(bin))
			if asymmErr != False :
				for process in asymmErr :
					file.write('\t%f' % histos_tmp[process].GetBinErrorUp(bin))
				for process in asymmErr :
					file.write('\t%f' % histos_tmp[process].GetBinErrorLow(bin))
			if bin_width :
				file.write('\t%f' % histos[processes[0]].GetBinWidth(bin))
			file.write('%s\n' % lumi_str)
		if last_bin :
			file.write('%f\t%f' % (histos[processes[0]].GetXaxis().GetBinUpEdge(nbins), histos[processes[0]].GetXaxis().GetBinCenter(nbins)))
			for process in processes :
				file.write('\t%f' % histos[process].GetBinContent(nbins))
			for process in processes :
				file.write('\t%f' % histos[process].GetBinError(nbins))
			if asymmErr != False :
				for process in asymmErr :
					file.write('\t%f' % histos_tmp[process].GetBinErrorUp(nbins))
				for process in asymmErr :
					file.write('\t%f' % histos_tmp[process].GetBinErrorLow(nbins))
			if bin_width :
				file.write('\t%f' % histos[processes[0]].GetBinWidth(nbins))
			file.write('%s\n' % lumi_str)

	# save histogram as CSV table
	if not asymmErr : asymmErr = []
	entries = {}
	for process in processes :
		# sanity check
		if (
				histos[process].GetNbinsX()                                          != histos[processes[0]].GetNbinsX() or
				histos[process].GetXaxis().GetBinLowEdge(1)                          != histos[processes[0]].GetXaxis().GetBinLowEdge(1) or
				histos[process].GetXaxis().GetBinUpEdge(histos[process].GetNbinsX()) != histos[processes[0]].GetXaxis().GetBinUpEdge(histos[processes[0]].GetNbinsX())
				) :
			print '[ERROR] histograms %s and %s have inconsistent binning!' % (process, processes[0])
			sys.exit(1)

		# get arrays from histogram
		entries.update(get_arraysFromHisto(histo = histos[process], var_str = var, data_str = process, last_bin = last_bin, x_range = x_range))
		if not 'stat' in entries.keys() :
			entries['stat']     = np.zeros(len(entries['binlow']))
			entries['stat_err'] = np.zeros(len(entries['binlow']))
			if lumi != '' :
				entries['stat'    ][0] = lumi
				entries['stat_err'][0] = 0.026 * lumi

	# save CSV table
	columns = ['binlow', var] + processes + ['%s_err' % process for process in processes] + ['%s_uerr' % process for process in asymmErr] + ['%s_derr' % process for process in asymmErr] + ['stat',]
	tables.write_CSVTable(entries, columns, filename, path)


def get_arraysFromHisto(histo, var_str = 'bin_centre', data_str = 'bin_content', last_bin = False, x_range = None) :
	histo_tmp = histo.Clone()
	histo_tmp.SetBinErrorOption(ROOT.TH1.kPoisson)
	entries = {}
	entries['binlow'            ] = []
	entries[             var_str] = []
	entries['%s'      % data_str] = []
	entries['%s_err'  % data_str] = []
	entries['%s_uerr' % data_str] = []
	entries['%s_derr' % data_str] = []
	if x_range == None :
		first_bin = 1
		nbins = histo.GetNbinsX()
	else :
		first_bin = histo.FindBin(x_range[0])
		nbins     = histo.FindBin(x_range[1])
	for ibin in range(first_bin, nbins+1) :
		entries['binlow'            ].append(histo.GetXaxis().GetBinLowEdge (ibin))
		entries[             var_str].append(histo.GetXaxis().GetBinCenter  (ibin))
		entries['%s'      % data_str].append(histo           .GetBinContent (ibin))
		entries['%s_err'  % data_str].append(histo           .GetBinError   (ibin))
		entries['%s_uerr' % data_str].append(histo_tmp       .GetBinErrorUp (ibin))
		entries['%s_derr' % data_str].append(histo_tmp       .GetBinErrorLow(ibin))
	if last_bin :
		entries['binlow'            ].append(histo.GetXaxis().GetBinUpEdge  (nbins))
		entries[             var_str].append(float('nan')                          )
		entries['%s'      % data_str].append(histo           .GetBinContent (nbins))
		entries['%s_err'  % data_str].append(histo           .GetBinError   (nbins))
		entries['%s_uerr' % data_str].append(histo_tmp       .GetBinErrorUp (nbins))
		entries['%s_derr' % data_str].append(histo_tmp       .GetBinErrorLow(nbins))
	for entry in entries :
		entries[entry] = np.array(entries[entry])
	return entries


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
