#! /usr/bin/python
import os, sys, commands, subprocess, math
import ROOT

ROOT.gROOT.ProcessLine( \
	"struct Entry{ \
		Double_t signif; \
	};"
	)
from ROOT import Entry

def getSignificance(rootfile):
	rootFile = ROOT.TFile.Open(rootfile);
	tree = rootFile.Get("limit");
	valueStruct = Entry()
	tree.SetBranchAddress("limit",ROOT.AddressOf(valueStruct,"signif"))
	tree.GetEntry(0)
	return valueStruct.signif

def observed_significance(datacard, opt):
	print 'combine -M ProfileLikelihood --signif '+datacard+' '+opt
	os.system('combine -M ProfileLikelihood --signif '+datacard+' '+opt)
	obs_signif = getSignificance("higgsCombineTest.ProfileLikelihood.mH120.root")
	print 'observed significance ', obs_signif
	return obs_signif

def expected_significance(datacard, opt):
	print 'combine -M ProfileLikelihood --significance '+datacard+' -t -1 --expectSignal=1'+' '+opt
	os.system('combine -M ProfileLikelihood --significance '+datacard+' -t -1 --expectSignal=1'+' '+opt)
	exp_signif = getSignificance("higgsCombineTest.ProfileLikelihood.mH120.root")
	print 'expected significance ', exp_signif
	return exp_signif

def signal_strength(datacard, opt):
	print 'combine -M MaxLikelihoodFit '+datacard+' '+opt
	os.system('combine -M MaxLikelihoodFit '+datacard+' '+opt)
	rootFile = ROOT.TFile.Open("mlfit.root")
	fit_s = rootFile.Get("fit_s")
#	fit_s.Print()
#	RooRealVar *rf = dynamic_cast<RooRealVar*>(res_s->floatParsFinal().find(r->GetName()));
	rf = fit_s.floatParsFinal().find("r")
	bestFitVal = rf.getVal()
	hiErr = +(rf.getMax("err68") - bestFitVal)
	loErr = -(rf.getMin("err68") - bestFitVal)
	maxError = max(max(hiErr, loErr), rf.getError())
	if (abs(hiErr) < 0.001*maxError):
		hiErr = -bestFitVal + rf.getMax()
	if (abs(loErr) < 0.001*maxError):
		loErr = +bestFitVal - rf.getMin()

	
	print 'Best fit r: ', rf.getVal(), '  ', -loErr, '/+', +hiErr, '  (68% CL)'
	
	return (rf.getVal(), loErr, hiErr)
	
#	print 'Best fit ', rf.getVal()

def printLaTeX(obsSignif, obsPValue, expSignif, expPValue, sigStrength, loStat, hiStat, loSyst, hiSyst, channel):
	xsec = 232.
	print '\\newcommand{\\ttWExpSignificance'+channel+'}{%.2f}' % expSignif
	print '\\newcommand{\\ttWExpPValue'+channel+'}      {%.2f}' % expPValue
	print '\\newcommand{\\ttWSignificance'+channel+'}   {%.2f}' % obsSignif
	print '\\newcommand{\\ttWPValue'+channel+'}         {%.2f}' % obsPValue
	print '\\newcommand{\\ttWCrossSection'+channel+'}   {%.0f  \\,\\,\\,^{+%.0f}_{-%.0f} \\,\\,\\,\\mathrm{(stat)} \\,\\,\\,  ^{+%.0f}_{-%.0f}\\,\\,\\, \\mathrm{(syst)} \\,\\,\\,  \\mathrm{fb}}' % (sigStrength*xsec, hiStat*xsec, loStat*xsec, hiSyst*xsec, loSyst*xsec)
#	print '\\renewcommand{\\ttWCrossSection}   {', sigStrength*xsec, '  \\,\\,\\,^{+', hiStat*xsec, '}_{-', loStat*xsec, '} \\,\\,\\,\\mathrm{(stat)} \\,\\,\\,  ^{+', hiSyst*xsec, '}_{-', loSyst*xsec, '}\\,\\,\\, \mathrm{(syst)} \\,\\,\\,  \\mathrm{pb}}'

def main(args):
	if ('--help' in args) or ('-h' in args):
		print 'add usage!'
		sys.exit()
	if ('-d' in args) or ('--datacard' in args):
		obsSignif = observed_significance(str(args[args.index('-d')+1]),"")
		expSignif = expected_significance(str(args[args.index('-d')+1]),"")
		obsPValue = observed_significance(str(args[args.index('-d')+1]),"--pvalue")
		expPValue = expected_significance(str(args[args.index('-d')+1]),"--pvalue")
		(sigStrength , loSystStat, hiSystStat) = signal_strength(str(args[args.index('-d')+1]),"")
		(sigStrength2, loStat    , hiStat    ) = signal_strength(str(args[args.index('-d')+1]),"-S 0")
#		(sigStrength2, loStat    , hiStat    ) = signal_strength(str(args[args.index('-d')+1]),"--justFit --profilingMode=none")
		loSyst = math.sqrt(loSystStat*loSystStat - loStat*loStat)
		hiSyst = math.sqrt(hiSystStat*hiSystStat - hiStat*hiStat)
		channel = ''
		if ('-c' in args) or ('--channel' in args):
			channel = str(args[args.index('-c')+1])
		printLaTeX(obsSignif, obsPValue, expSignif, expPValue, sigStrength, loStat, hiStat, loSyst, hiSyst, channel)


print '[status] starting...'
main(sys.argv)
print '[status] ...done'


