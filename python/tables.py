#! /usr/bin/python
import helper
import time
import ttvplot
import math


def make_ObsPredTable(path, results) :
#		table_name = 'datacard_ssdl_ttW_' + results['Normal'][charge][chan].chan_str + '.txt'
	table_name = 'ObsPredTable.tex'
	table_path = path + 'IntPredictions/'
	helper.mkdir(table_path)
	pl = ttvplot.ttvplot(table_path, '2L', TeX_switch = True)
	print '[status] writing %s' % table_name
	with open(table_path + table_name, 'w') as file :
		timestamp = time.asctime()
		file.write('%!TEX root = ../AN-12-445.tex\n')
		file.write('%=========================================================================================\n')
		file.write('% Observation and predictions table for ttW analysis, same-sign channel, subchannels\n')
		file.write('%% Generated on: %s\n' % str(timestamp))
		file.write('%-----------------------------------------------------------------------------------------\n')
		file.write('\n\n')
		file.write('\\providecommand{\\ttbar}{$t\\bar{t}$}\n')
		file.write('\\providecommand{\\ttw}{\\ttbar{}W}\n')
		file.write('\\providecommand{\\ttz}{\\ttbar{}Z}\n')
		file.write('\\providecommand{\\Pep}{e^+}\n')
		file.write('\\providecommand{\\Pem}{e^-}\n')
		file.write('\\providecommand{\\Pgmp}{\\mu^+}\n')
		file.write('\\providecommand{\\Pgmm}{\\mu^-}\n')
		file.write('\n\n')
	##	fOUTSTREAM << "%% Format is m+m+, e+m+, e+e+, m-m-, e-m-, e-e-" << endl;
		file.write('\\begin{tabular}{l|r@{$\\,\\pm\\,$}l|r@{$\\,\\pm\\,$}l|r@{$\\,\\pm\\,$}l|r@{$\\,\\pm\\,$}l|r@{$\\,\\pm\\,$}l|r@{$\\,\\pm\\,$}l}\n\\hline \\hline\n')
		file.write('                  &\\multicolumn{2}{c|}{$\\Pgmp\\Pgmp$}&\\multicolumn{2}{c|}{$\\Pep\\Pgmp$}&\\multicolumn{2}{c|}{$\\Pep\\Pep$} &\\multicolumn{2}{c|}{$\\Pgmm\\Pgmm$}&\\multicolumn{2}{c|}{$\\Pem\\Pgmm$}&\\multicolumn{2}{c}{$\\Pem\\Pem$} \\\\\n\\hline\n')
		file.write('%18s & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f \\\\\n' % (
			pl.get_processName('fake'),
			results['++']['mm'].fake, results['++']['mm'].fake_err,
			results['++']['em'].fake, results['++']['em'].fake_err,
			results['++']['ee'].fake, results['++']['ee'].fake_err,
			results['--']['mm'].fake, results['--']['mm'].fake_err,
			results['--']['em'].fake, results['--']['em'].fake_err,
			results['--']['ee'].fake, results['--']['ee'].fake_err))
	##//	fOUTSTREAM3 << Form("Double Fakes   & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \n",
	##//					   nff_mm, sqrt(FR->getMMNffEStat()*FR->getMMNffEStat()+nff_mm*nff_mm*FakeESyst2),
	##//					   nff_em, sqrt(FR->getEMNffEStat()*FR->getEMNffEStat()+nff_em*nff_em*FakeESyst2),
	##//					   nff_ee, sqrt(FR->getEENffEStat()*FR->getEENffEStat()+nff_ee*nff_ee*FakeESyst2),
	##//					   nff_em + nff_mm + nff_ee, sqrt(FR->getTotDoubleEStat()*FR->getTotDoubleEStat() + nDF*nDF*FakeESyst2));
	##//	fOUTSTREAM3 << Form("Single Fakes   & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \n",
	##//					   npf_mm,          sqrt(FR->getMMNpfEStat()   *FR->getMMNpfEStat()    +  npf_mm*npf_mm*FakeESyst2),
	##//					   npf_em + nfp_em, sqrt(FR->getEMSingleEStat()*FR->getEMSingleEStat() + (npf_em+nfp_em)*(npf_em+nfp_em)*FakeESyst2),
	##//					   npf_ee,          sqrt(FR->getEENpfEStat()   *FR->getEENpfEStat()    +  npf_ee*npf_ee*FakeESyst2),
	##//					   npf_em + nfp_em + npf_mm + npf_ee, sqrt(FR->getTotSingleEStat()*FR->getTotSingleEStat() + nSF*nSF*FakeESyst2));
		file.write('%18s & \\multicolumn{2}{c|}{-}        & %13.1f & %13.1f & %13.1f & %13.1f & \\multicolumn{2}{c|}{-}        & %13.1f & %13.1f & %13.1f & %13.1f \\\\\n' % (
			pl.get_processName('chmid'),
			results['++']['em'].cmid, results['++']['em'].cmid_err,
			results['++']['ee'].cmid, results['++']['ee'].cmid_err,
			results['--']['em'].cmid, results['--']['em'].cmid_err,
			results['--']['ee'].cmid, results['--']['ee'].cmid_err))
		file.write('%18s & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f \\\\\n' % (
			pl.get_processName('rare'),
			results['++']['mm'].rare, results['++']['mm'].rare_err,
			results['++']['em'].rare, results['++']['em'].rare_err,
			results['++']['ee'].rare, results['++']['ee'].rare_err,
			results['--']['mm'].rare, results['--']['mm'].rare_err,
			results['--']['em'].rare, results['--']['em'].rare_err,
			results['--']['ee'].rare, results['--']['ee'].rare_err))
		file.write('%18s & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f \\\\\n' % (
			pl.get_processName('wz'),
			results['++']['mm'].wz, results['++']['mm'].wz_err,
			results['++']['em'].wz, results['++']['em'].wz_err,
			results['++']['ee'].wz, results['++']['ee'].wz_err,
			results['--']['mm'].wz, results['--']['mm'].wz_err,
			results['--']['em'].wz, results['--']['em'].wz_err,
			results['--']['ee'].wz, results['--']['ee'].wz_err))
		file.write('%18s & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f \\\\\n' % (
			pl.get_processName('ttz'),
			results['++']['mm'].ttz, results['++']['mm'].ttz_err,
			results['++']['em'].ttz, results['++']['em'].ttz_err,
			results['++']['ee'].ttz, results['++']['ee'].ttz_err,
			results['--']['mm'].ttz, results['--']['mm'].ttz_err,
			results['--']['em'].ttz, results['--']['em'].ttz_err,
			results['--']['ee'].ttz, results['--']['ee'].ttz_err))
		file.write('\\hline\n')
	##//	if (separateTTH) {
	##//		fOUTSTREAM3 << Form("ttH Prod.      & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \\hline \n",
	##//							tth_nt2_mm, sqrt(tth_nt2_mm_e1 + RareESyst2*tth_nt2_mm*tth_nt2_mm),
	##//							tth_nt2_em, sqrt(tth_nt2_em_e1 + RareESyst2*tth_nt2_em*tth_nt2_em),
	##//							tth_nt2_ee, sqrt(tth_nt2_ee_e1 + RareESyst2*tth_nt2_ee*tth_nt2_ee),
	##//							tth_nt2_ee + tth_nt2_mm + tth_nt2_em, sqrt(tth_nt2_mm_e1 + tth_nt2_ee_e1 + tth_nt2_em_e1 + RareESyst2*(tth_nt2_ee + tth_nt2_mm + tth_nt2_em)*(tth_nt2_ee + tth_nt2_mm + tth_nt2_em)));
	##//	}
		file.write('%18s & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f \\\\\n' % (
			'Total Bkg.',
			results['++']['mm'].tot, results['++']['mm'].tot_err,
			results['++']['em'].tot, results['++']['em'].tot_err,
			results['++']['ee'].tot, results['++']['ee'].tot_err,
			results['--']['mm'].tot, results['--']['mm'].tot_err,
			results['--']['em'].tot, results['--']['em'].tot_err,
			results['--']['ee'].tot, results['--']['ee'].tot_err))
		file.write('{\\bf %12s} &  \\multicolumn{2}{c|}{\\bf %3d} &  \\multicolumn{2}{c|}{\\bf %3d} &  \\multicolumn{2}{c|}{\\bf %3d} &  \\multicolumn{2}{c|}{\\bf %3d} &  \\multicolumn{2}{c|}{\\bf %3d} &  \\multicolumn{2}{c}{\\bf %3d}  \\\\\n' % (
			pl.get_processName('obs'),
			results['++']['mm'].obs,
			results['++']['em'].obs,
			results['++']['ee'].obs,
			results['--']['mm'].obs,
			results['--']['em'].obs,
			results['--']['ee'].obs))
		file.write('%18s & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f \\\\\n' % (
			'Obs. $-$ Tot. Bkg.',
			results['++']['mm'].obs - results['++']['mm'].tot, results['++']['mm'].tot_err,
			results['++']['em'].obs - results['++']['em'].tot, results['++']['em'].tot_err,
			results['++']['ee'].obs - results['++']['ee'].tot, results['++']['ee'].tot_err,
			results['--']['mm'].obs - results['--']['mm'].tot, results['--']['mm'].tot_err,
			results['--']['em'].obs - results['--']['em'].tot, results['--']['em'].tot_err,
			results['--']['ee'].obs - results['--']['ee'].tot, results['--']['ee'].tot_err))
		file.write('\\hline\n')
#		file.write('%18s & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f \\\\\n' % (
#			pl.get_processName('ttw') + ' (exp.)',
#			results['++']['mm'].ttw, results['++']['mm'].ttw_err,
#			results['++']['em'].ttw, results['++']['em'].ttw_err,
#			results['++']['ee'].ttw, results['++']['ee'].ttw_err,
#			results['--']['mm'].ttw, results['--']['mm'].ttw_err,
#			results['--']['em'].ttw, results['--']['em'].ttw_err,
#			results['--']['ee'].ttw, results['--']['ee'].ttw_err))
		file.write('%18s & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f & %13.1f \\\\\n' % (
			pl.get_processName('ttw') + ' (exp.)',
			results['++']['mm'].ttw, math.sqrt(results['++']['mm'].ttw_staterr*results['++']['mm'].ttw_staterr + 0.08*0.08*results['++']['mm'].ttw*results['++']['mm'].ttw),
			results['++']['em'].ttw, math.sqrt(results['++']['em'].ttw_staterr*results['++']['em'].ttw_staterr + 0.08*0.08*results['++']['em'].ttw*results['++']['em'].ttw),
			results['++']['ee'].ttw, math.sqrt(results['++']['ee'].ttw_staterr*results['++']['ee'].ttw_staterr + 0.08*0.08*results['++']['ee'].ttw*results['++']['ee'].ttw),
			results['--']['mm'].ttw, math.sqrt(results['--']['mm'].ttw_staterr*results['--']['mm'].ttw_staterr + 0.08*0.08*results['--']['mm'].ttw*results['--']['mm'].ttw),
			results['--']['em'].ttw, math.sqrt(results['--']['em'].ttw_staterr*results['--']['em'].ttw_staterr + 0.08*0.08*results['--']['em'].ttw*results['--']['em'].ttw),
			results['--']['ee'].ttw, math.sqrt(results['--']['ee'].ttw_staterr*results['--']['ee'].ttw_staterr + 0.08*0.08*results['--']['ee'].ttw*results['--']['ee'].ttw)))
		file.write('\\hline \\hline\n')
		file.write('\\end{tabular}\n')
		file.write('\n\n')


def make_SystTable(path, results, chan, charge) :
	'''
	takes a nested dictionary of result objects as input:
	results[SYSTFLAG][CHARGE][FLAVOR]

	writes systematics table
	'''

	channel = 'Int'
	table_name = 'SystTable'+channel+'.tex'
	table_path = path + 'IntPredictions/'
	helper.mkdir(table_path)
	print '[status] writing %s' % table_name
	with open(table_path + table_name, 'w') as file :
		timestamp = time.asctime()
		file.write('%!TEX root = ../AN-12-445.tex\n')
		file.write('%=========================================================================================\n')
		file.write('% Systematics table for ttW analysis, same-sign channel, subchannels\n')
		file.write('%% Generated on: %s\n' % str(timestamp))
		file.write('%-----------------------------------------------------------------------------------------\n')
		file.write('\n\n')
		file.write('\\begin{tabular}{l|r@{$\\,/\\,$}l}\n')
		file.write('\\hline\\hline\n')
		file.write('Systematic & \\multicolumn{2}{c}{$\\Delta\\epsilon$ [\\%]} \\\\\n')
		file.write('\\hline\n')

		if 'JetUp' in results and 'JetDown' in results :
			file.write('JES up/down            & %6.1f & %6.1f \\\\\n' % (
				100. * (results['JetUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['JetDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.)
				))
		else :
			print '[WARNING] JetUp/Down systematic not found!'

		if 'JetSmearUp' in results and 'JetSmearDown' in results :
			file.write('JER up/down            & %6.1f & %6.1f \\\\\n' % (
				100. * (results['JetSmearUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['JetSmearDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.)
				))
		else :
			print '[WARNING] JetSmearUp/Down systematic not found!'

		if 'BUp' in results and 'BDown' in results :
			file.write('b-tag up/down          & %6.1f & %6.1f \\\\\n' % (
				100. * (results['BUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['BDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.)
				))
		else :
			print '[WARNING] BUp/Down systematic not found!'

		if 'LepUp' in results and 'LepDown' in results :
			file.write('Lepton scale up/down   & %6.1f & %6.1f \\\\\n' % (
				100. * (results['LepUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['LepDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.)
				))
		else :
			print '[WARNING] LepUp/Down systematic not found!'

		if 'MuUp' in results and 'MuDown' in results :
			file.write('Muon scale up/down     & %6.1f & %6.1f \\\\\n' % (
				100. * (results['MuUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['MuDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.)
				))
		else :
			print '[WARNING] MuUp/Down systematic not found!'

		if 'ElUp' in results and 'ElDown' in results :
			file.write('Electron scale up/down & %6.1f & %6.1f \\\\\n' % (
				100. * (results['ElUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['ElDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.)
				))
		else :
			print '[WARNING] ElUp/Down systematic not found!'

		if 'PileupUp' in results and 'PileupDown' in results :
			file.write('Pileup up/down         & %6.1f & %6.1f \\\\\n' % (
				100. * (results['PileupUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['PileupDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.)
				))
		else :
			print '[WARNING] PileupUp/Down systematic not found!'

		file.write('\\hline\\hline')
		file.write('\\end{tabular}')
