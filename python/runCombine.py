#! /usr/bin/python
import os, sys, commands, subprocess, math
import ROOT
import helper


def getSignificance(filepath) :
	file = ROOT.TFile.Open(filepath)
	tree = file.Get('limit')
	tree.GetEntry(0)
	return tree.limit


def observed_significance(datacard, opt = '') :
	print '[status] calculating observed significance..'
	combineCommand = 'combine -M ProfileLikelihood --signif %s %s' % (datacard, opt)
	print combineCommand
	os.system(combineCommand)
	obs_signif = getSignificance("higgsCombineTest.ProfileLikelihood.mH120.root")
	print 'observed significance: %5.2f' % obs_signif
	return obs_signif


def expected_significance(datacard, opt = '') :
	print '[status] calculating expected significance..'
	combineCommand = 'combine -M ProfileLikelihood --significance %s -t -1 --expectSignal=1 %s' % (datacard, opt)
	print combineCommand
	os.system(combineCommand)
	exp_signif = getSignificance("higgsCombineTest.ProfileLikelihood.mH120.root")
	print 'expected significance: %5.2f' % exp_signif
	return exp_signif


def signal_strength(datacard, opt = '') :
	print '[status] calculating signal strength..'
	combineCommand = 'combine -M MaxLikelihoodFit %s %s' % (datacard, opt)
	print combineCommand
	os.system(combineCommand)
	file = ROOT.TFile.Open('mlfit.root')
	fit_s = file.Get('fit_s')
#	fit_s.Print()
#	RooRealVar *rf = dynamic_cast<RooRealVar*>(res_s->floatParsFinal().find(r->GetName()));
	rf = fit_s.floatParsFinal().find('r')
	bestFitVal = rf.getVal()
	hiErr = +(rf.getMax('err68') - bestFitVal)
	loErr = -(rf.getMin('err68') - bestFitVal)
	maxError = max(max(hiErr, loErr), rf.getError())
	if (abs(hiErr) < 0.001*maxError) :
		hiErr = -bestFitVal + rf.getMax()
	if (abs(loErr) < 0.001*maxError) :
		loErr = +bestFitVal - rf.getMin()
	print 'Best fit r: %6.3f -%5.3f/+%5.3f (68%% CL)' % (rf.getVal(), loErr, hiErr)
	return (rf.getVal(), loErr, hiErr)


def printLaTeX(obsSignif, obsPValue, expSignif, expPValue, sigStrength, loStat, hiStat, loSyst, hiSyst, channel) :
	xsec = 206.
	print '\\newcommand{\\ttWExpSignificance'+channel+'}{%.2f}' % expSignif
	print '\\newcommand{\\ttWExpPValue'+channel+'}      {%.2f}' % expPValue
	print '\\newcommand{\\ttWSignificance'+channel+'}   {%.2f}' % obsSignif
	print '\\newcommand{\\ttWPValue'+channel+'}         {%.2f}' % obsPValue
	print '\\newcommand{\\ttWCrossSection'+channel+'}   {%.0f  \\,\\,\\,^{+%.0f}_{-%.0f} \\,\\,\\,\\mathrm{(stat)} \\,\\,\\,  ^{+%.0f}_{-%.0f}\\,\\,\\, \\mathrm{(syst)} \\,\\,\\,  \\mathrm{fb}}' % (sigStrength*xsec, hiStat*xsec, loStat*xsec, hiSyst*xsec, loSyst*xsec)
#	print '\\renewcommand{\\ttWCrossSection}   {', sigStrength*xsec, '  \\,\\,\\,^{+', hiStat*xsec, '}_{-', loStat*xsec, '} \\,\\,\\,\\mathrm{(stat)} \\,\\,\\,  ^{+', hiSyst*xsec, '}_{-', loSyst*xsec, '}\\,\\,\\, \mathrm{(syst)} \\,\\,\\,  \\mathrm{pb}}'
	print 'sigma = %.0f +%.0f -%.0f (stat) +%.0f -%.0f (syst)' % (sigStrength*xsec, hiStat*xsec, loStat*xsec, hiSyst*xsec, loSyst*xsec)
	print 'sigma = %.1f +%.1f -%.1f (tot)'                     % (sigStrength*xsec, math.sqrt(hiStat**2 + hiSyst**2)*xsec, math.sqrt(loStat**2 + loSyst**2)*xsec)
	print '\nrounded:'
	print 'sigma = %s +%s -%s (stat) +%s -%s (syst)' % helper.get_roundedNumberErrors(sigStrength*xsec, [hiStat*xsec, loStat*xsec, hiSyst*xsec, loSyst*xsec])
	print '{%s}{+%s}{-%s}{+%s}{-%s}' % helper.get_roundedNumberErrors(sigStrength*xsec, [hiStat*xsec, loStat*xsec, hiSyst*xsec, loSyst*xsec])
#	print '\\renewcommand{\\ttWCrossSection}   {', sigStrength*xsec, '  \\,\\,\\,^{+', hiStat*xsec, '}_{-', loStat*xsec, '} \\,\\,\\,\\mathrm{(stat)} \\,\\,\\,  ^{+', hiSyst*xsec, '}_{-', loSyst*xsec, '}\\,\\,\\, \mathrm{(syst)} \\,\\,\\,  \\mathrm{pb}}'


def main(args) :
	if ('--help' in args) or ('-h' in args) :
		print 'add usage!'
		sys.exit()
	if ('-d' in args) or ('--datacard' in args) :
		obsSignif = observed_significance(str(args[args.index('-d')+1]), '')
		expSignif = expected_significance(str(args[args.index('-d')+1]), '')
		obsPValue = observed_significance(str(args[args.index('-d')+1]), '--pvalue')
		expPValue = expected_significance(str(args[args.index('-d')+1]), '--pvalue')
		(sigStrength , loSystStat, hiSystStat) = signal_strength(str(args[args.index('-d')+1]), '')
		(sigStrength2, loStat    , hiStat    ) = signal_strength(str(args[args.index('-d')+1]), '-S 0')
#		(sigStrength2, loStat    , hiStat    ) = signal_strength(str(args[args.index('-d')+1]), '--justFit --profilingMode=none')
		loSyst = math.sqrt(loSystStat*loSystStat - loStat*loStat)
		hiSyst = math.sqrt(hiSystStat*hiSystStat - hiStat*hiStat)
		channel = ''
		if ('-c' in args) or ('--channel' in args):
			channel = str(args[args.index('-c')+1])
		printLaTeX(obsSignif, obsPValue, expSignif, expPValue, sigStrength, loStat, hiStat, loSyst, hiSyst, channel)


if __name__ == '__main__' :
	print '[status] starting..'
	main(sys.argv)
	print '[status] done'
