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
			'Obs.\\ $-$ Tot.\\ Bkg.',
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


def make_SystTable(path, results, chan, charge, systematics) :
	'''
	takes a nested dictionary of result objects as input:
	results[SYSTFLAG][CHARGE][FLAVOR]

	writes systematics table
	'''

	channel = 'Int'
	table_name = 'SystTable'+channel+'.tex'
	table_path = path + 'IntPredictions/'
	helper.mkdir(table_path)
	pl = ttvplot.ttvplot(table_path, '2L', TeX_switch = True, short_names = True)
	print '[status] writing %s' % table_name
	with open(table_path + table_name, 'w') as file :
		timestamp = time.asctime()
		file.write('%!TEX root = ../AN-12-445.tex\n')
		file.write('%=========================================================================================\n')
		file.write('% Systematics table for ttW analysis, same-sign channel, subchannels\n')
		file.write('%% Generated on: %s\n' % str(timestamp))
		file.write('%-----------------------------------------------------------------------------------------\n')
		file.write('\n\n')
		file.write('\\providecommand{\\ttbar}{$t\\bar{t}$}\n')
		file.write('\\providecommand{\\ttw}{\\ttbar{}W}\n')
		file.write('\\providecommand{\\ttz}{\\ttbar{}Z}\n')
		file.write('\n\n')
		file.write('\\begin{tabular}{l|c|ccccc}\n')
		file.write('\\hline\\hline\n')
		file.write('{\\bf Source} & {\\bf Signal} [\\%] & \\multicolumn{5}{c}{{\\bf Backgrounds} [\\%]} \\\\\n')
		file.write('& $\\Delta\\varepsilon_s$ & %s & %s & %s & %s & %s \\\\\n' % (
			pl.get_processName('ttz'  ),
			pl.get_processName('fake' ),
			pl.get_processName('chmid'),
			pl.get_processName('wz'   ),
			pl.get_processName('rare' )
			))
		file.write('\\hline\n')

		file.write('%18s & - & %6.1f & %6.1f & %6.1f & %6.1f & %6.1f \\\\\n' % (
			'Bkg.\\ Pred.',
			100. * results['Normal'][charge][chan].ttz_err  / results['Normal'][charge][chan].ttz ,
			100. * results['Normal'][charge][chan].fake_err / results['Normal'][charge][chan].fake,
			100. * results['Normal'][charge][chan].cmid_err / results['Normal'][charge][chan].cmid,
			100. * results['Normal'][charge][chan].wz_err   / results['Normal'][charge][chan].wz  ,
			100. * results['Normal'][charge][chan].rare_err / results['Normal'][charge][chan].rare
			))

		if 'lumi' in systematics :
			file.write('%18s & %6.1f & %6.1f & - & - & %6.1f & %6.1f \\\\\n' % (
				'Luminosity',
				100. * (systematics['lumi']-1.),
				100. * (systematics['lumi']-1.),
				100. * (systematics['lumi']-1.),
				100. * (systematics['lumi']-1.)
				))
		else :
			print '[WARNING] Luminosity systematic not found!'

		if 'JetUp' in results and 'JetDown' in results :
			file.write('%18s & %6.1f / %6.1f & %6.1f / %6.1f & - & - & %6.1f / %6.1f & %6.1f / %6.1f \\\\\n' % (
				'JES',
				100. * (results['JetUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['JetDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['JetUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['JetDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['JetUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['JetDown'][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['JetUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare - 1.),
				100. * (results['JetDown'][charge][chan].rare / results['Normal'][charge][chan].rare - 1.)
				))
		else :
			print '[WARNING] JetUp/Down systematic not found!'

		if 'JetSmearUp' in results and 'JetSmearDown' in results :
			file.write('%18s & %6.1f / %6.1f & %6.1f / %6.1f & - & - & %6.1f / %6.1f & %6.1f / %6.1f \\\\\n' % (
				'JER',
				100. * (results['JetSmearUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['JetSmearDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['JetSmearUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['JetSmearDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['JetSmearUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['JetSmearDown'][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['JetSmearUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare - 1.),
				100. * (results['JetSmearDown'][charge][chan].rare / results['Normal'][charge][chan].rare - 1.)
				))
		else :
			print '[WARNING] JetSmearUp/Down systematic not found!'

		if 'BUp' in results and 'BDown' in results :
			file.write('%18s & %6.1f / %6.1f & %6.1f / %6.1f & - & - & %6.1f / %6.1f & %6.1f / %6.1f \\\\\n' % (
				'b-tagging',
				100. * (results['BUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['BDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['BUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['BDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['BUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['BDown'][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['BUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare - 1.),
				100. * (results['BDown'][charge][chan].rare / results['Normal'][charge][chan].rare - 1.)
				))
		else :
			print '[WARNING] BUp/Down systematic not found!'

		if 'ltrig' in systematics :
			file.write('%18s & %6.1f & %6.1f & - & - & %6.1f & %6.1f \\\\\n' % (
				'Lept.\\ Trig.',
				100. * (systematics['ltrig']-1.),
				100. * (systematics['ltrig']-1.),
				100. * (systematics['ltrig']-1.),
				100. * (systematics['ltrig']-1.)
				))
		else :
			print '[WARNING] Trigger systematic not found!'

		if 'LepUp' in results and 'LepDown' in results :
			file.write('%18s & %6.1f / %6.1f & %6.1f / %6.1f & - & - & %6.1f / %6.1f & %6.1f / %6.1f \\\\\n' % (
				'Lepton SF',
				100. * (results['LepUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['LepDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['LepUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['LepDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['LepUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['LepDown'][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['LepUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare - 1.),
				100. * (results['LepDown'][charge][chan].rare / results['Normal'][charge][chan].rare - 1.)
				))
		else :
			print '[WARNING] LepUp/Down systematic not found!'

		if 'MuUp' in results and 'MuDown' in results :
			file.write('%18s & %6.1f / %6.1f & %6.1f / %6.1f & - & - & %6.1f / %6.1f & %6.1f / %6.1f \\\\\n' % (
				'Muon SF',
				100. * (results['MuUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['MuDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['MuUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['MuDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['MuUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['MuDown'][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['MuUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare - 1.),
				100. * (results['MuDown'][charge][chan].rare / results['Normal'][charge][chan].rare - 1.)
				))
		else :
			print '[WARNING] MuUp/Down systematic not found!'

		if 'ElUp' in results and 'ElDown' in results :
			file.write('%18s & %6.1f / %6.1f & %6.1f / %6.1f & - & - & %6.1f / %6.1f & %6.1f / %6.1f \\\\\n' % (
				'Electron SF',
				100. * (results['ElUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['ElDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['ElUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['ElDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['ElUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['ElDown'][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['ElUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare - 1.),
				100. * (results['ElDown'][charge][chan].rare / results['Normal'][charge][chan].rare - 1.)
				))
		else :
			print '[WARNING] ElUp/Down systematic not found!'

		if 'PileupUp' in results and 'PileupDown' in results :
			file.write('%18s & %6.1f / %6.1f & %6.1f / %6.1f & - & - & %6.1f / %6.1f & %6.1f / %6.1f \\\\\n' % (
				'Pileup',
				100. * (results['PileupUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['PileupDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['PileupUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['PileupDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz - 1.),
				100. * (results['PileupUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['PileupDown'][charge][chan].wz   / results['Normal'][charge][chan].wz   - 1.),
				100. * (results['PileupUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare - 1.),
				100. * (results['PileupDown'][charge][chan].rare / results['Normal'][charge][chan].rare - 1.)
				))
		else :
			print '[WARNING] PileupUp/Down systematic not found!'

		if 'scale_up' in systematics and 'scale_dn' in systematics :
			file.write('%18s & %6.1f / %6.1f & %6.1f / %6.1f & - & - & %6.1f / %6.1f & %6.1f / %6.1f \\\\\n' % (
				'$Q^2$',
				100. * (systematics['scale_up']-1.),
				100. * (systematics['scale_dn']-1.),
				100. * (systematics['scale_up']-1.),
				100. * (systematics['scale_dn']-1.),
				100. * (systematics['scale_up']-1.),
				100. * (systematics['scale_dn']-1.),
				100. * (systematics['scale_up']-1.),
				100. * (systematics['scale_dn']-1.)
				))
		else :
			print '[WARNING] Scale systematic not found!'

		if 'tmass_up' in systematics and 'tmass_dn' in systematics :
			file.write('%18s & %6.1f / %6.1f & %6.1f / %6.1f & - & - & - & - \\\\\n' % (
				'$m_{\\text{top}}$',
				100. * (systematics['tmass_up']-1.),
				100. * (systematics['tmass_dn']-1.),
				100. * (systematics['tmass_up']-1.),
				100. * (systematics['tmass_dn']-1.)
				))
		else :
			print '[WARNING] Top mass systematic not found!'

		if 'gen' in systematics :
			file.write('%18s & %6.1f & - & - & - & - & - \\\\\n' % (
				'Generator',
				100. * (systematics['gen']-1.)
				))
		else :
			print '[WARNING] Generator systematic not found!'

		if 'pdf' in systematics and 'wz_pdf' in systematics :
			file.write('%18s & %6.1f & - & - & - & %6.1f & - \\\\\n' % (
				'PDF',
				100. * (systematics['pdf']-1.),
				100. * (systematics['wz_pdf']-1.)
				))
		else :
			print '[WARNING] PDF systematic not found!'

		file.write('\\hline\\hline\n')
		file.write('\\end{tabular}')


def make_YieldsTable(res, systematics) :
	'''print all observations and predictions'''

	# PRINTOUT
	print "-------------------------------------------------------------------------------------------------------------------------------"
	print "                 |               Mu/Mu               |                E/Mu               |                E/E                ||"
	print "         YIELDS  |   Ntt  |   Ntl  |   Nlt  |   Nll  |   Ntt  |   Ntl  |   Nlt  |   Nll  |   Ntt  |   Ntl  |   Nlt  |   Nll  ||"
	print "-------------------------------------------------------------------------------------------------------------------------------"
	print "%16s & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f" % ("Data",
		res['mm'].nt2 ,
		res['mm'].nt10,
		res['mm'].nt01,
		res['mm'].nt0 ,
		res['em'].nt2 ,
		res['em'].nt10,
		res['em'].nt01,
		res['em'].nt0 ,
		res['ee'].nt2 ,
		res['ee'].nt10,
		res['ee'].nt01,
		res['ee'].nt0 )



	print "  Fake Predictions:"
	print "------------------------------------------------------------------------------------------------------"
	print "                 |            Mu/Mu          |           El/Mu           |            El/El          |"
	print "------------------------------------------------------------------------------------------------------"
	print " Npp             |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
		res['mm'].npp, res['mm'].npp_staterr, res['mm'].npp_systerr,
		res['em'].npp, res['em'].npp_staterr, res['em'].npp_systerr,
		res['ee'].npp, res['ee'].npp_staterr, res['ee'].npp_systerr)
	print " Npf             |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
		res['mm'].npf, res['mm'].npf_staterr, res['mm'].npf_systerr,
		res['em'].npf, res['em'].npf_staterr, res['em'].npf_systerr,
		res['ee'].npf, res['ee'].npf_staterr, res['ee'].npf_systerr)
	print " Nfp             |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
		res['mm'].nfp, res['mm'].nfp_staterr, res['mm'].nfp_systerr,
		res['em'].nfp, res['em'].nfp_staterr, res['em'].nfp_systerr,
		res['ee'].nfp, res['ee'].nfp_staterr, res['ee'].nfp_systerr)
	print " Nff             |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
		res['mm'].nff, res['mm'].nff_staterr, res['mm'].nff_systerr,
		res['em'].nff, res['em'].nff_staterr, res['em'].nff_systerr,
		res['ee'].nff, res['ee'].nff_staterr, res['ee'].nff_systerr)
	print "------------------------------------------------------------------------------------------------------"
	print " Total Fakes     |", "%5.1f +/- %5.1f           | %5.1f +/- %5.1f           | %5.1f +/- %5.1f           |" % (
		res['mm'].fake, res['mm'].fake_err,
		res['em'].fake, res['em'].fake_err,
		res['ee'].fake, res['ee'].fake_err)
	print "------------------------------------------------------------------------------------------------------"
	print " (Value +/- E_stat +/- E_syst) "
	print "//////////////////////////////////////////////////////////////////////////////////////////"
	print " ChMisID         |", "                          | %5.1f +/- %5.1f           | %5.1f +/- %5.1f           |" % (
		res['em'].cmid, res['em'].cmid_err,
		res['ee'].cmid, res['ee'].cmid_err)
	print "------------------------------------------------------------------------------------------------------"

	print 'rares:'

	for s in res['al'].rares :
		print "%16s || %5.2f +/- %5.2f +/- %5.2f || %5.2f +/- %5.2f +/- %5.2f || %5.2f +/- %5.2f +/- %5.2f || %5.2f +/- %5.2f +/- %5.2f ||" % (s, 
			res['mm'].rares[s], res['mm'].rares_staterr[s], systematics['rare']*res['mm'].rares[s],
			res['em'].rares[s], res['em'].rares_staterr[s], systematics['rare']*res['em'].rares[s],
			res['ee'].rares[s], res['ee'].rares_staterr[s], systematics['rare']*res['ee'].rares[s],
			res['al'].rares[s], res['al'].rares_staterr[s], systematics['rare']*res['al'].rares[s])
	print "----------------------------------------------------------------------------------------------"
	print "%16s || %5.2f +/- %5.2f +/- %5.2f || %5.2f +/- %5.2f +/- %5.2f || %5.2f +/- %5.2f +/- %5.2f || %5.2f +/- %5.2f +/- %5.2f ||" % ('Total', 
		res['mm'].rare, res['mm'].rare_staterr, systematics['rare']*res['mm'].rare,
		res['em'].rare, res['em'].rare_staterr, systematics['rare']*res['em'].rare,
		res['ee'].rare, res['ee'].rare_staterr, systematics['rare']*res['ee'].rare,
		res['al'].rare, res['al'].rare_staterr, systematics['rare']*res['al'].rare)


#		print "----------------------------------------------------------------------------------------------"
#		print "       SUMMARY   ||         Mu/Mu         ||         E/Mu          ||          E/E          ||"
#		print "=============================================================================================="
#		print "%16s || %5.2f +/- %5.2f +/- %5.2f || %5.2f +/- %5.2f +/- %5.2f || %5.2f +/- %5.2f +/- %5.2f ||\n" % ("pred. fakes",
#			self.nF_mm, FR.getMMTotEStat(), self.FakeESyst*nF_mm,
#			self.nF_em, FR.getEMTotEStat(), self.FakeESyst*nF_em,
#			self.nF_ee, FR.getEETotEStat(), self.FakeESyst*nF_ee)
#		print "%16s ||                       || %5.2f +/- %5.2f +/- %5.2f || %5.2f +/- %5.2f +/- %5.2f ||\n" % ("pred. chmisid",
#			nt2_em_chmid, nt2_em_chmid_e1, nt2_em_chmid_e2, nt2_ee_chmid, nt2_ee_chmid_e1, nt2_ee_chmid_e2)
#
	print "----------------------------------------------------------------------------------------------"


def make_CutsTable(path, selections) :
	'''
	writes cuts table

	selections is a list of tuples: (eff, cuts)
	'''

	table_name = 'cuts.tex'
	table_path = path
	if not table_path.endswith('/') : table_path += '/'
	helper.mkdir(table_path)
#	pl = ttvplot.ttvplot(table_path, '2L', TeX_switch = True, short_names = True)
	print '[status] writing %s' % table_name
	ncuts = len(selections[0][1])
	with open(table_path + table_name, 'w') as file :
		timestamp = time.asctime()
		file.write('%!TEX root = ../AN-12-445.tex\n')
		file.write('%=========================================================================================\n')
		file.write('% Cuts table for ttW analysis, same-sign channel, subchannels\n')
		file.write('%% Generated on: %s\n' % str(timestamp))
		file.write('%-----------------------------------------------------------------------------------------\n')
		file.write('\n\n')
		file.write('\\begin{tabular}{l')
		for i in range(ncuts) :
			file.write('|cc')
		file.write('}\n')
		file.write('\\hline\\hline\n')
		file.write('{\\bf Eff.}')
		for var in selections[0][1] :
#			file.write(' & \\multicolumn{2}{|c}{\\bf %s}' % pl.get_varName(var))
			file.write(' & \\multicolumn{2}{|c}{\\bf %s}' % var)
		file.write(' \\\\\n')
		for i in range(ncuts) :
			file.write(' & min & max')
		file.write(' \\\\\n')
		file.write('\\hline\n')
		for (eff, cuts) in selections :
			file.write('%3.0f \\%%' % eff)
			for var in cuts :
				file.write(' & %5.0f & %5.0e' % (cuts[var][0], cuts[var][1]))
			file.write(' \\\\\n')
		file.write('\\hline\\hline\n')
		file.write('\\end{tabular}')