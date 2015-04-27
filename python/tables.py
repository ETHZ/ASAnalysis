#! /usr/bin/python
import helper
import time
import ttvStyle
import math


def make_ObsPredTable(path, results) :
#		table_name = 'datacard_ssdl_ttW_' + results['Normal'][charge][chan].chan_str + '.txt'
	table_name = 'ObsPredTable.tex'
	table_path = path + 'IntPredictions/'
	rnd = None
	float_digits = None
	helper.mkdir(table_path)
	pl = ttvStyle.ttvStyle(TeX_switch = True)
	print '[status] writing %s' % table_name
	with open(table_path + table_name, 'w') as file :
		timestamp = time.asctime()
		file.write('%!TEX root = ../../Dissertation.tex\n')
		file.write('%=========================================================================================\n')
		file.write('% Observation and predictions table for ttW analysis, same-sign channel, subchannels\n')
		file.write('%% Generated on: %s\n' % str(timestamp))
		file.write('%-----------------------------------------------------------------------------------------\n')
		file.write('\n\n')
		file.write('\\providecommand{\\ttbar}{\\ensuremath{t\\bar{t}}}\n')
		file.write('\\providecommand{\\ttw}  {\\ensuremath{\\ttbar{}W}}\n')
		file.write('\\providecommand{\\ttz}  {\\ensuremath{\\ttbar{}Z}}\n')
		file.write('\\providecommand{\\wz}   {\\ensuremath{WZ}}\n')
		file.write('\\providecommand{\\Pep}  {\\ensuremath{e^+}}\n')
		file.write('\\providecommand{\\Pem}  {\\ensuremath{e^-}}\n')
		file.write('\\providecommand{\\PGmp} {\\ensuremath{\\mu^+}}\n')
		file.write('\\providecommand{\\PGmm} {\\ensuremath{\\mu^-}}\n')
		file.write('\n\n')
		file.write('\\sisetup{\n')
		file.write('\tseparate-uncertainty\n')
		file.write('}\n')
		file.write('\n\n')
		file.write('\\begin{tabular}{\n')
		file.write('\tl\n')
		file.write('\tS[table-number-alignment = center, table-format = 1.2, table-figures-uncertainty = 2]\n')
		file.write('\tS[table-number-alignment = center, table-format = 1.2, table-figures-uncertainty = 2]\n')
		file.write('\tS[table-number-alignment = center, table-format = 1.2, table-figures-uncertainty = 2]\n')
		file.write('\tS[table-number-alignment = center, table-format = 1.2, table-figures-uncertainty = 2]\n')
		file.write('\tS[table-number-alignment = center, table-format = 1.2, table-figures-uncertainty = 2]\n')
		file.write('\tS[table-number-alignment = center, table-format = 1.2, table-figures-uncertainty = 2]\n')
		file.write('}\n')
		file.write('\t\\toprule\n')
		file.write('\t                     &  {\\PGmp\\PGmp} &  {\\Pep\\PGmp}  &   {\\Pep\\Pep}  &  {\\PGmm\\PGmm} &  {\\Pem\\PGmm}  &   {\\Pem\\Pem}  \\\\\n')
		file.write('\t\\midrule\n')
#		file.write('\t%-20s & %5.1f +- %3.1f & %5.1f +- %3.1f & %5.1f +- %3.1f & %5.1f +- %3.1f & %5.1f +- %3.1f & %5.1f +- %3.1f \\\\\n' % (
		file.write('\t%-20s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s \\\\\n' % (
			(pl.get_processName('fake'),) +
			helper.get_roundedNumber(results['++']['mm'].fake, results['++']['mm'].fake_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['em'].fake, results['++']['em'].fake_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['ee'].fake, results['++']['ee'].fake_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['mm'].fake, results['--']['mm'].fake_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['em'].fake, results['--']['em'].fake_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['ee'].fake, results['--']['ee'].fake_err, rnd, float_digits)))
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
		file.write('\t%-20s &      {--}     & %5s +- %4s & %5s +- %4s &      {--}     & %5s +- %4s & %5s +- %4s \\\\\n' % (
			(pl.get_processName('chmid'),) +
			helper.get_roundedNumber(results['++']['em'].cmid, results['++']['em'].cmid_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['ee'].cmid, results['++']['ee'].cmid_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['em'].cmid, results['--']['em'].cmid_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['ee'].cmid, results['--']['ee'].cmid_err, rnd, float_digits)))
		file.write('\t%-20s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s \\\\\n' % (
			(pl.get_processName('rare'),) +
			helper.get_roundedNumber(results['++']['mm'].rare, results['++']['mm'].rare_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['em'].rare, results['++']['em'].rare_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['ee'].rare, results['++']['ee'].rare_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['mm'].rare, results['--']['mm'].rare_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['em'].rare, results['--']['em'].rare_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['ee'].rare, results['--']['ee'].rare_err, rnd, float_digits)))
		file.write('\t%-20s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s \\\\\n' % (
			(pl.get_processName('wz'),) +
			helper.get_roundedNumber(results['++']['mm'].wz, results['++']['mm'].wz_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['em'].wz, results['++']['em'].wz_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['ee'].wz, results['++']['ee'].wz_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['mm'].wz, results['--']['mm'].wz_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['em'].wz, results['--']['em'].wz_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['ee'].wz, results['--']['ee'].wz_err, rnd, float_digits)))
		file.write('\t%-20s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s \\\\\n' % (
			(pl.get_processName('ttz'),) +
			helper.get_roundedNumber(results['++']['mm'].ttz, results['++']['mm'].ttz_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['em'].ttz, results['++']['em'].ttz_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['ee'].ttz, results['++']['ee'].ttz_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['mm'].ttz, results['--']['mm'].ttz_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['em'].ttz, results['--']['em'].ttz_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['ee'].ttz, results['--']['ee'].ttz_err, rnd, float_digits)))
	##//	if (separateTTH) {
	##//		fOUTSTREAM3 << Form("ttH Prod.      & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f & %5.1f &$\\pm$ %5.1f \\\\ \\midrule \n",
	##//							tth_nt2_mm, sqrt(tth_nt2_mm_e1 + RareESyst2*tth_nt2_mm*tth_nt2_mm),
	##//							tth_nt2_em, sqrt(tth_nt2_em_e1 + RareESyst2*tth_nt2_em*tth_nt2_em),
	##//							tth_nt2_ee, sqrt(tth_nt2_ee_e1 + RareESyst2*tth_nt2_ee*tth_nt2_ee),
	##//							tth_nt2_ee + tth_nt2_mm + tth_nt2_em, sqrt(tth_nt2_mm_e1 + tth_nt2_ee_e1 + tth_nt2_em_e1 + RareESyst2*(tth_nt2_ee + tth_nt2_mm + tth_nt2_em)*(tth_nt2_ee + tth_nt2_mm + tth_nt2_em)));
	##//	}
		file.write('\t%-20s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s \\\\\n' % (
			('Total background',) +
			helper.get_roundedNumber(results['++']['mm'].tot, results['++']['mm'].tot_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['em'].tot, results['++']['em'].tot_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['ee'].tot, results['++']['ee'].tot_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['mm'].tot, results['--']['mm'].tot_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['em'].tot, results['--']['em'].tot_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['ee'].tot, results['--']['ee'].tot_err, rnd, float_digits)))
#		file.write('\t%-20s & %5.1f +- %3.1f & %5.1f +- %3.1f & %5.1f +- %3.1f & %5.1f +- %3.1f & %5.1f +- %3.1f & %5.1f +- %3.1f \\\\\n' % (
#			'Obs.\\ $-$ Tot.\\ Bkg.',
#			results['++']['mm'].obs - results['++']['mm'].tot, results['++']['mm'].tot_err,
#			results['++']['em'].obs - results['++']['em'].tot, results['++']['em'].tot_err,
#			results['++']['ee'].obs - results['++']['ee'].tot, results['++']['ee'].tot_err,
#			results['--']['mm'].obs - results['--']['mm'].tot, results['--']['mm'].tot_err,
#			results['--']['em'].obs - results['--']['em'].tot, results['--']['em'].tot_err,
#			results['--']['ee'].obs - results['--']['ee'].tot, results['--']['ee'].tot_err))
		file.write('\t\\midrule\n')
#		file.write('\t%-20s & %5.1f +- %3.1f & %5.1f +- %3.1f & %5.1f +- %3.1f & %5.1f +- %3.1f & %5.1f +- %3.1f & %5.1f +- %3.1f \\\\\n' % (
#			pl.get_processName('ttw') + ' (exp.)',
#			results['++']['mm'].ttw, results['++']['mm'].ttw_err,
#			results['++']['em'].ttw, results['++']['em'].ttw_err,
#			results['++']['ee'].ttw, results['++']['ee'].ttw_err,
#			results['--']['mm'].ttw, results['--']['mm'].ttw_err,
#			results['--']['em'].ttw, results['--']['em'].ttw_err,
#			results['--']['ee'].ttw, results['--']['ee'].ttw_err))
		file.write('\t%-20s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s \\\\\n' % (
			(pl.get_processName('ttw') + '{} (expected)',) +
			helper.get_roundedNumber(results['++']['mm'].ttw, math.sqrt(results['++']['mm'].ttw_staterr*results['++']['mm'].ttw_staterr + 0.08*0.08*results['++']['mm'].ttw*results['++']['mm'].ttw), rnd, float_digits) +
			helper.get_roundedNumber(results['++']['em'].ttw, math.sqrt(results['++']['em'].ttw_staterr*results['++']['em'].ttw_staterr + 0.08*0.08*results['++']['em'].ttw*results['++']['em'].ttw), rnd, float_digits) +
			helper.get_roundedNumber(results['++']['ee'].ttw, math.sqrt(results['++']['ee'].ttw_staterr*results['++']['ee'].ttw_staterr + 0.08*0.08*results['++']['ee'].ttw*results['++']['ee'].ttw), rnd, float_digits) +
			helper.get_roundedNumber(results['--']['mm'].ttw, math.sqrt(results['--']['mm'].ttw_staterr*results['--']['mm'].ttw_staterr + 0.08*0.08*results['--']['mm'].ttw*results['--']['mm'].ttw), rnd, float_digits) +
			helper.get_roundedNumber(results['--']['em'].ttw, math.sqrt(results['--']['em'].ttw_staterr*results['--']['em'].ttw_staterr + 0.08*0.08*results['--']['em'].ttw*results['--']['em'].ttw), rnd, float_digits) +
			helper.get_roundedNumber(results['--']['ee'].ttw, math.sqrt(results['--']['ee'].ttw_staterr*results['--']['ee'].ttw_staterr + 0.08*0.08*results['--']['ee'].ttw*results['--']['ee'].ttw), rnd, float_digits)))
		file.write('\t%-20s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s & %5s +- %4s \\\\\n' % (
			('Total expected',) +
			helper.get_roundedNumber(results['++']['mm'].tot_exp, results['++']['mm'].tot_exp_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['em'].tot_exp, results['++']['em'].tot_exp_err, rnd, float_digits) +
			helper.get_roundedNumber(results['++']['ee'].tot_exp, results['++']['ee'].tot_exp_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['mm'].tot_exp, results['--']['mm'].tot_exp_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['em'].tot_exp, results['--']['em'].tot_exp_err, rnd, float_digits) +
			helper.get_roundedNumber(results['--']['ee'].tot_exp, results['--']['ee'].tot_exp_err, rnd, float_digits)))
		file.write('\t{\\bf %-14s} &  {\\bf %6d} &  {\\bf %6d} &  {\\bf %6d} &  {\\bf %6d} &  {\\bf %6d} &  {\\bf %6d} \\\\\n' % (
			pl.get_processName('obs'),
			results['++']['mm'].obs,
			results['++']['em'].obs,
			results['++']['ee'].obs,
			results['--']['mm'].obs,
			results['--']['em'].obs,
			results['--']['ee'].obs))
		file.write('\t\\bottomrule\n')
		file.write('\\end{tabular}\n')


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
	pl = ttvStyle.ttvStyle(TeX_switch = True, short_names = True)
	print '[status] writing %s' % table_name
	with open(table_path + table_name, 'w') as file :
		timestamp = time.asctime()
		file.write('%!TEX root = ../../Dissertation.tex\n')
		file.write('%=========================================================================================\n')
		file.write('% Systematics table for ttW analysis, same-sign channel, subchannels\n')
		file.write('%% Generated on: %s\n' % str(timestamp))
		file.write('%-----------------------------------------------------------------------------------------\n')
		file.write('\n\n')
		file.write('\\providecommand{\\ttbar}{$t\\bar{t}$}\n')
		file.write('\\providecommand{\\ttw}{\\ttbar{}W}\n')
		file.write('\\providecommand{\\ttz}{\\ttbar{}Z}\n')
		file.write('\\providecommand{\\wz}{WZ}\n')
		file.write('\n\n')
		file.write('\\begin{tabular}{\n')
		file.write('\tl\n')
		file.write('\tS[table-format = -1.1, retain-explicit-plus]@{ / }\n')
		file.write('\tS[table-format = -1.1, retain-explicit-plus]\n')
		file.write('\tS[table-format = -1.1, retain-explicit-plus]@{ / }\n')
		file.write('\tS[table-format = -1.1, retain-explicit-plus]\n')
		file.write('\tS[table-format = 2.1]\n')
		file.write('\tS[table-format = 2.1]\n')
		file.write('\tS[table-format = -1.1, retain-explicit-plus]@{ / }\n')
		file.write('\tS[table-format = -1.1, retain-explicit-plus]\n')
		file.write('\tS[table-format = -1.1, retain-explicit-plus]@{ / }\n')
		file.write('\tS[table-format = -1.1, retain-explicit-plus]\n')
		file.write('}\n')
		file.write('\t\\toprule\n')
		file.write('\t{\\bf Source}       & \\multicolumn{2}{c}{{\\bf Signal} (\\si{\\percent})} & \\multicolumn{8}{c}{{\\bf Backgrounds} (\\si{\\percent})} \\\\\n')
		file.write('\t\\cmidrule(lr){2-3}\n')
		file.write('\t\\cmidrule(lr){4-11}\n')
		file.write('\t                   & \\multicolumn{2}{c}{$\\Delta\\varepsilon_s$} & \\multicolumn{2}{c}{%s} & {%s} & {%s} & \\multicolumn{2}{c}{%s} & \\multicolumn{2}{c}{%s} \\\\\n' % (
			pl.get_processName('ttz'  ),
			pl.get_processName('fake' ),
			pl.get_processName('chmid'),
			pl.get_processName('wz'   ),
			pl.get_processName('rare' )
			))
		file.write('\t\\midrule\n')

		file.write('\t%-18s & \\multicolumn{2}{c}{-}               & \\multicolumn{2}{c}{\\tablenum{%4.1f}} & %4.1f & %4.1f & \\multicolumn{2}{c}{\\tablenum{%4.1f}} & \\multicolumn{2}{c}{\\tablenum{%4.1f}} \\\\\n' % (
			'Bkg.\\ Pred.',
			100. * results['Normal'][charge][chan].ttz_err  / results['Normal'][charge][chan].ttz ,
			100. * results['Normal'][charge][chan].fake_err / results['Normal'][charge][chan].fake,
			100. * results['Normal'][charge][chan].cmid_err / results['Normal'][charge][chan].cmid,
			100. * results['Normal'][charge][chan].wz_err   / results['Normal'][charge][chan].wz  ,
			100. * results['Normal'][charge][chan].rare_err / results['Normal'][charge][chan].rare
			))

		if 'lumi' in systematics :
			file.write('\t%-18s & \\multicolumn{2}{c}{\\tablenum{%4.1f}} & \\multicolumn{2}{c}{\\tablenum{%4.1f}} &  {-} &  {-} & \\multicolumn{2}{c}{\\tablenum{%4.1f}} & \\multicolumn{2}{c}{\\tablenum{%4.1f}} \\\\\n' % (
				'Luminosity',
				100. * (systematics['lumi']-1.),
				100. * (systematics['lumi']-1.),
				100. * (systematics['lumi']-1.),
				100. * (systematics['lumi']-1.)
				))
		else :
			print '[WARNING] Luminosity systematic not found!'

		if 'JetUp' in results and 'JetDown' in results :
			file.write('\t%-18s &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            &  {-} &  {-} &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            \\\\\n' % (
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
			file.write('\t%-18s &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            &  {-} &  {-} &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            \\\\\n' % (
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
			file.write('\t%-18s &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            &  {-} &  {-} &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            \\\\\n' % (
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
			file.write('\t%-18s & \\multicolumn{2}{c}{\\tablenum{%4.1f}} & \\multicolumn{2}{c}{\\tablenum{%4.1f}} &  {-} &  {-} & \\multicolumn{2}{c}{\\tablenum{%4.1f}} & \\multicolumn{2}{c}{\\tablenum{%4.1f}} \\\\\n' % (
				'Lept.\\ Trig.',
				100. * (systematics['ltrig']-1.),
				100. * (systematics['ltrig']-1.),
				100. * (systematics['ltrig']-1.),
				100. * (systematics['ltrig']-1.)
				))
		else :
			print '[WARNING] Trigger systematic not found!'

		if 'LepUp' in results and 'LepDown' in results :
			file.write('\t%-18s &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            &  {-} &  {-} &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            \\\\\n' % (
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
			file.write('\t%-18s &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            &  {-} &  {-} &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            \\\\\n' % (
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
			file.write('\t%-18s &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            &  {-} &  {-} &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            \\\\\n' % (
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
			file.write('\t%-18s &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            &  {-} &  {-} &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            \\\\\n' % (
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
			file.write('\t%-18s &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            &  {-} &  {-} &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            \\\\\n' % (
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
			file.write('\t%-18s &            %+5.1f & %+5.1f            &            %+5.1f & %+5.1f            &  {-} &  {-} & \\multicolumn{2}{c}{-}               & \\multicolumn{2}{c}{-}               \\\\\n' % (
				'$m_{\\text{top}}$',
				100. * (systematics['tmass_up']-1.),
				100. * (systematics['tmass_dn']-1.),
				100. * (systematics['tmass_up']-1.),
				100. * (systematics['tmass_dn']-1.)
				))
		else :
			print '[WARNING] Top mass systematic not found!'

		if 'gen' in systematics :
			file.write('\t%-18s & \\multicolumn{2}{c}{\\tablenum{%4.1f}} & \\multicolumn{2}{c}{-}               &  {-} &  {-} & \\multicolumn{2}{c}{-}               & \\multicolumn{2}{c}{-}               \\\\\n' % (
				'Generator',
				100. * (systematics['gen']-1.)
				))
		else :
			print '[WARNING] Generator systematic not found!'

		if 'pdf' in systematics and 'wz_pdf' in systematics :
			file.write('\t%-18s & \\multicolumn{2}{c}{\\tablenum{%4.1f}} & \\multicolumn{2}{c}{-}               &  {-} &  {-} & \\multicolumn{2}{c}{\\tablenum{%4.1f}} & \\multicolumn{2}{c}{-}               \\\\\n' % (
				'PDF',
				100. * (systematics['pdf']-1.),
				100. * (systematics['wz_pdf']-1.)
				))
		else :
			print '[WARNING] PDF systematic not found!'

		file.write('\t\\bottomrule\n')
		file.write('\\end{tabular}\n')


def make_YieldsTable(path, res, systematics, suffix = '') :
	'''print all observations and predictions'''

	if suffix != '' and not suffix.startswith('_') :
		suffix = '_' + suffix

	# data yields
	table_name = 'DataYieldsTable%s.tex' % suffix
	table_path = path + 'IntPredictions/'
	helper.mkdir(table_path)
	pl = ttvStyle.ttvStyle(TeX_switch = True)
	print '[status] writing %s' % table_name
	with open(table_path + table_name, 'w') as file :
		timestamp = time.asctime()
		file.write('%!TEX root = ../../Dissertation.tex\n')
		file.write('\n')
		file.write(providecommands())
		file.write('\n\n')
		file.write('\\begin{tabular}{\n')
		file.write('\tl\n')
		file.write('\tS[table-number-alignment = center, table-format = 2.0]\n')
		file.write('\tS[table-number-alignment = center, table-format = 2.0]\n')
		file.write('\tS[table-number-alignment = center, table-format = 2.0]\n')
		file.write('}\n')
		file.write('\t\\toprule\n')
		file.write('\tYields & {\\PGm\\PGm} & {\\Pe\\PGm} & {\\Pe\\Pe} \\\\\n')
		file.write('\t\\midrule\n')
		file.write('\t%-6s &   %8.0f &  %8.0f & %8.0f \\\\\n' % (
			'\\ntt',
			res['mm'].nt2,
			res['em'].nt2,
			res['ee'].nt2))
		file.write('\t%-6s &   %8.0f &  %8.0f & %8.0f \\\\\n' % (
			'\\ntl',
			res['mm'].nt10,
			res['em'].nt10,
			res['ee'].nt10))
		file.write('\t%-6s &   %8.0f &  %8.0f & %8.0f \\\\\n' % (
			'\\nlt',
			res['mm'].nt01,
			res['em'].nt01,
			res['ee'].nt01))
		file.write('\t%-6s &   %8.0f &  %8.0f & %8.0f \\\\\n' % (
			'\\nll',
			res['mm'].nt0,
			res['em'].nt0,
			res['ee'].nt0))
		file.write('\t\\bottomrule\n')
		file.write('\\end{tabular}\n')

	# predictions
	table_name = 'PredTable%s.tex' % suffix
	table_path = path + 'IntPredictions/'
	helper.mkdir(table_path)
	pl = ttvStyle.ttvStyle(TeX_switch = True)
	print '[status] writing %s' % table_name
	with open(table_path + table_name, 'w') as file :
		timestamp = time.asctime()
		file.write('%!TEX root = ../../Dissertation.tex\n')
		file.write('\n')
		file.write(providecommands())
		file.write('\n\n')
		file.write('\\sisetup{separate-uncertainty}\n')
		file.write('\n\n')
		file.write('\\begin{tabular}{\n')
		file.write('\tl\n')
		file.write('\tS[table-number-alignment = center, table-format = 2.1, table-figures-uncertainty = 2]\n')
		file.write('\tS[table-number-alignment = center, table-format = 2.1, table-figures-uncertainty = 2]\n')
		file.write('\tS[table-number-alignment = center, table-format = 2.1, table-figures-uncertainty = 2]\n')
		file.write('}\n')
		file.write('\t\\toprule\n')
		file.write('\tPrediction   &     {\\PGm\\PGm}     &      {\\Pe\\PGm}     &      {\\Pe\\Pe}      \\\\\n')
		# fakes
		file.write('\t\\midrule\n')
		file.write('\t%-12s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('\\npp',) +
			helper.get_roundedNumber(res['mm'].npp, res['mm'].npp_staterr) + #res['mm'].npp_systerr,
			helper.get_roundedNumber(res['em'].npp, res['em'].npp_staterr) + #res['em'].npp_systerr,
			helper.get_roundedNumber(res['ee'].npp, res['ee'].npp_staterr)))#, res['ee'].npp_systerr))
		file.write('\t%-12s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('\\npf',) +
			helper.get_roundedNumber(res['mm'].npf, res['mm'].npf_staterr) + #res['mm'].npf_systerr,
			helper.get_roundedNumber(res['em'].npf, res['em'].npf_staterr) + #res['em'].npf_systerr,
			helper.get_roundedNumber(res['ee'].npf, res['ee'].npf_staterr)))#, res['ee'].npf_systerr))
		file.write('\t%-12s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('\\nfp',) +
			helper.get_roundedNumber(res['mm'].nfp, res['mm'].nfp_staterr) + #res['mm'].nfp_systerr,
			helper.get_roundedNumber(res['em'].nfp, res['em'].nfp_staterr) + #res['em'].nfp_systerr,
			helper.get_roundedNumber(res['ee'].nfp, res['ee'].nfp_staterr)))#, res['ee'].nfp_systerr))
		file.write('\t%-12s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('\\nff',) +
			helper.get_roundedNumber(res['mm'].nff, res['mm'].nff_staterr) + #res['mm'].nff_systerr,
			helper.get_roundedNumber(res['em'].nff, res['em'].nff_staterr) + #res['em'].nff_systerr,
			helper.get_roundedNumber(res['ee'].nff, res['ee'].nff_staterr)))#, res['ee'].nff_systerr))
		file.write('\t\\midrule\n')
		file.write('\t%-12s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('Total Fakes',) +
			helper.get_roundedNumber(res['mm'].fake, res['mm'].fake_staterr) +
			helper.get_roundedNumber(res['em'].fake, res['em'].fake_staterr) +
			helper.get_roundedNumber(res['ee'].fake, res['ee'].fake_staterr)))
		# charge mis-ID
		file.write('\t\\midrule\n')
		file.write('\t%-12s &                    & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('Charge MisID',) +
			helper.get_roundedNumber(res['em'].cmid, res['em'].cmid_staterr) +
			helper.get_roundedNumber(res['ee'].cmid, res['ee'].cmid_staterr)))
		# rares
		file.write('\t\\midrule\n')
		for s in res['al'].rares :
			file.write('\t%-12s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
				(s,) +
				helper.get_roundedNumber(res['mm'].rares[s], res['mm'].rares_staterr[s]) +# systematics['rare']*res['mm'].rares[s],
				helper.get_roundedNumber(res['em'].rares[s], res['em'].rares_staterr[s]) +# systematics['rare']*res['em'].rares[s],
				helper.get_roundedNumber(res['ee'].rares[s], res['ee'].rares_staterr[s])))#, systematics['rare']*res['ee'].rares[s],
		file.write('\t\\midrule\n')
		file.write('\t%-12s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('Rares',) +
			helper.get_roundedNumber(res['mm'].rare, res['mm'].rare_staterr) +# systematics['rare']*res['mm'].rare,
			helper.get_roundedNumber(res['em'].rare, res['em'].rare_staterr) +# systematics['rare']*res['em'].rare,
			helper.get_roundedNumber(res['ee'].rare, res['ee'].rare_staterr)))#, systematics['rare']*res['ee'].rare,
		file.write('\t\\midrule\n')
		file.write('\t%-12s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('\\wz',) +
			helper.get_roundedNumber(res['mm'].wz, res['mm'].wz_staterr) +
			helper.get_roundedNumber(res['em'].wz, res['em'].wz_staterr) +
			helper.get_roundedNumber(res['ee'].wz, res['ee'].wz_staterr)))
		file.write('\t\\midrule\n')
		file.write('\t%-12s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('\\ttz',) +
			helper.get_roundedNumber(res['mm'].ttz, res['mm'].ttz_staterr) +
			helper.get_roundedNumber(res['em'].ttz, res['em'].ttz_staterr) +
			helper.get_roundedNumber(res['ee'].ttz, res['ee'].ttz_staterr)))
		print "res['em'].ttz", res['em'].ttz
		file.write('\t\\bottomrule\n')
		file.write('\\end{tabular}\n')


def make_MCYieldsTable(path, res, suffix = '') :
	'''print all MC yields'''

	if suffix != '' and not suffix.startswith('_') :
		suffix = '_' + suffix

	# data yields
	table_name = 'MCYieldsTable%s.tex' % suffix
	table_path = path + 'IntPredictions/'
	helper.mkdir(table_path)
	pl = ttvStyle.ttvStyle(TeX_switch = True)
	print '[status] writing %s' % table_name
	with open(table_path + table_name, 'w') as file :
		timestamp = time.asctime()
		file.write('%!TEX root = ../../Dissertation.tex\n')
		file.write('\n')
		file.write(providecommands())
		file.write('\n\n')
		file.write('\\begin{tabular}{\n')
		file.write('\tl\n')
		file.write('\tS[table-number-alignment = center, table-format = 2.1, table-figures-uncertainty = 2]\n')
		file.write('\tS[table-number-alignment = center, table-format = 2.1, table-figures-uncertainty = 2]\n')
		file.write('\tS[table-number-alignment = center, table-format = 2.1, table-figures-uncertainty = 2]\n')
		file.write('}\n')
		file.write('\t\\toprule\n')
		file.write('\tYields & {\\PGm\\PGm} & {\\Pe\\PGm} & {\\Pe\\Pe} \\\\\n')
		file.write('\t\\midrule\n')
		for s in res['al'].mc :
			file.write('\t%-14s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
				(s,) +
				helper.get_roundedNumber(res['mm'].mc[s], res['mm'].mc_staterr[s]) +# systematics['rare']*res['mm'].mc[s],
				helper.get_roundedNumber(res['em'].mc[s], res['em'].mc_staterr[s]) +# systematics['rare']*res['em'].mc[s],
				helper.get_roundedNumber(res['ee'].mc[s], res['ee'].mc_staterr[s])))#, systematics['rare']*res['ee'].mc[s],
		file.write('\t\\midrule\n')
		file.write('\t\\midrule\n')
		file.write('\t%-14s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('ttbar',) +
			helper.get_roundedNumber(res['mm'].ttbar, res['mm'].ttbar_staterr) +# systematics['rare']*res['mm'].mc[s],
			helper.get_roundedNumber(res['em'].ttbar, res['em'].ttbar_staterr) +# systematics['rare']*res['em'].mc[s],
			helper.get_roundedNumber(res['ee'].ttbar, res['ee'].ttbar_staterr)))#, systematics['rare']*res['ee'].mc[s],
		file.write('\t%-14s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('Single Top',) +
			helper.get_roundedNumber(res['mm'].singletop, res['mm'].singletop_staterr) +# systematics['rare']*res['mm'].mc[s],
			helper.get_roundedNumber(res['em'].singletop, res['em'].singletop_staterr) +# systematics['rare']*res['em'].mc[s],
			helper.get_roundedNumber(res['ee'].singletop, res['ee'].singletop_staterr)))#, systematics['rare']*res['ee'].mc[s],
		file.write('\t%-14s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('\\wjets',) +
			helper.get_roundedNumber(res['mm'].wjets, res['mm'].wjets_staterr) +# systematics['rare']*res['mm'].mc[s],
			helper.get_roundedNumber(res['em'].wjets, res['em'].wjets_staterr) +# systematics['rare']*res['em'].mc[s],
			helper.get_roundedNumber(res['ee'].wjets, res['ee'].wjets_staterr)))#, systematics['rare']*res['ee'].mc[s],
		file.write('\t%-14s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('\\zjets',) +
			helper.get_roundedNumber(res['mm'].zjets, res['mm'].zjets_staterr) +# systematics['rare']*res['mm'].mc[s],
			helper.get_roundedNumber(res['em'].zjets, res['em'].zjets_staterr) +# systematics['rare']*res['em'].mc[s],
			helper.get_roundedNumber(res['ee'].zjets, res['ee'].zjets_staterr)))#, systematics['rare']*res['ee'].mc[s],
		file.write('\t%-14s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('Irreducible',) +
			helper.get_roundedNumber(res['mm'].rare, res['mm'].rare_staterr) +# systematics['rare']*res['mm'].mc[s],
			helper.get_roundedNumber(res['em'].rare, res['em'].rare_staterr) +# systematics['rare']*res['em'].mc[s],
			helper.get_roundedNumber(res['ee'].rare, res['ee'].rare_staterr)))#, systematics['rare']*res['ee'].mc[s],
		file.write('\t%-14s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('\\wz',) +
			helper.get_roundedNumber(res['mm'].wz, res['mm'].wz_staterr) +# systematics['wz']*res['mm'].mc[s],
			helper.get_roundedNumber(res['em'].wz, res['em'].wz_staterr) +# systematics['wz']*res['em'].mc[s],
			helper.get_roundedNumber(res['ee'].wz, res['ee'].wz_staterr)))#, systematics['wz']*res['ee'].mc[s],
		file.write('\t%-14s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('\\ttz',) +
			helper.get_roundedNumber(res['mm'].ttz, res['mm'].ttz_staterr) +# systematics['\\ttz']*res['mm'].mc[s],
			helper.get_roundedNumber(res['em'].ttz, res['em'].ttz_staterr) +# systematics['\\ttz']*res['em'].mc[s],
			helper.get_roundedNumber(res['ee'].ttz, res['ee'].ttz_staterr)))#, systematics['\\ttz']*res['ee'].mc[s],
		file.write('\t\\midrule\n')
		file.write('\t%-14s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('\\ttw',) +
			helper.get_roundedNumber(res['mm'].ttw, res['mm'].ttw_staterr) +# systematics['ttw']*res['mm'].mc[s],
			helper.get_roundedNumber(res['em'].ttw, res['em'].ttw_staterr) +# systematics['ttw']*res['em'].mc[s],
			helper.get_roundedNumber(res['ee'].ttw, res['ee'].ttw_staterr)))#, systematics['ttw']*res['ee'].mc[s],
		file.write('\t%-14s & %8s +- %6s & %8s +- %6s & %8s +- %6s \\\\\n' % (
			('Others',) +
			helper.get_roundedNumber(res['mm'].other, res['mm'].other_staterr) +# systematics['others']*res['mm'].mc[s],
			helper.get_roundedNumber(res['em'].other, res['em'].other_staterr) +# systematics['others']*res['em'].mc[s],
			helper.get_roundedNumber(res['ee'].other, res['ee'].other_staterr)))#, systematics['others']*res['ee'].mc[s],
		file.write('\t\\bottomrule\n')
		file.write('\\end{tabular}\n')


def make_CutsTable(path, selections) :
	'''
	writes cuts table

	selections is a list of tuples: (eff, cuts)
	'''

	table_name = 'cuts.tex'
	table_path = path
	if not table_path.endswith('/') : table_path += '/'
	helper.mkdir(table_path)
#	pl = ttvStyle.ttvStyle(TeX_switch = True, short_names = True)
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


def make_FullOptTable(path, FoM, table, charge_str) :
	'''
	writes table with optimization results

	selections is a list of tuples: (eff, cuts)
	'''

	table_name = 'OptCutsTable_%s.tex' % charge_str
	table_path = path
	if not table_path.endswith('/') : table_path += '/'
	helper.mkdir(table_path)
	print '[status] writing %s' % table_name
	with open(table_path + table_name, 'w') as file :
		timestamp = time.asctime()
		file.write('%!TEX root = ../../Dissertation.tex\n')
		file.write('\n\n')
		file.write('\\begin{tabular}{')
		file.write('l') # column for ieff
		for cut in table[0]['cuts'] :
			file.write('|SS') # columns for cuts
		file.write('|S|') # column for efficiency
		for FoM in table[0]['results'] :
			for chan in table[0]['results'][FoM] :
				file.write('S') # columns for results
		file.write('}\n')
		file.write('\\hline\\hline\n')

		# first title row
		file.write('{\\bf Eff.}')
		for var in table[0]['cuts'] :
#			file.write(' & \\multicolumn{2}{c|}{\\bf %s}' % pl.get_varName(var))
			file.write(' & \\multicolumn{2}{c|}{\\bf %s}' % var)
		file.write(' & {$\\varepsilon$}')
		for FoM in table[0]['results'] :
			file.write(' & \\multicolumn{%d}{c|}{%s}' % (len(table[0]['results'][FoM]), FoM))
		file.write(' \\\\\n')

		# second title row
		for cut in table[0]['cuts'] :
			file.write('& {min} & {max} ')
		file.write('&')
		for FoM in table[0]['results'] :
			for chan in table[0]['results'][FoM] :
				file.write(' & {%s}' % chan)
		file.write(' \\\\\n')
		file.write('\\hline\n')

		# data rows
		for row in table :
			file.write('%3d' % row['ieff'])
			for var in row['cuts'] :
				file.write(' & %5.0f & %5.0e' % (row['cuts'][var][0], row['cuts'][var][1]))
			file.write(' & %4.1f +- %3.1f' % (row['eff'][0]*100., row['eff'][1]*100.))
			for FoM in row['results'] :
				for chan in row['results'][FoM] :
					file.write(' & %5.1f' % row['results'][FoM][chan][1])
			file.write(' \\\\\n')
		file.write('\\hline\\hline\n')
		file.write('\\end{tabular}')


def make_OptTable(path, FoM, table, charge_str) :
	'''
	writes table with optimization results

	selections is a list of tuples: (eff, cuts)
	'''

	table_name = 'OptCutsTable_%s.tex' % charge_str
	table_path = path
	if not table_path.endswith('/') : table_path += '/'
	helper.mkdir(table_path)
	print '[status] writing %s' % table_name
	with open(table_path + table_name, 'w') as file :
		timestamp = time.asctime()
		file.write('%!TEX root = ../../Dissertation.tex\n')
		file.write('\n\n')
		file.write('\\sisetup{\n')
		file.write('\tseparate-uncertainty,\n')
		file.write('\ttable-number-alignment = center,\n')
		file.write('\ttable-figures-integer  = 1,\n')
		file.write('\ttable-figures-decimal  = 0\n')
		file.write('}\n')
		file.write('\\begin{tabular}{')
#		file.write('l') # column for ieff
		for cut in table[0]['cuts'] :
			option = ''
			if table[0]['cuts'][cut][0] >= 1. :
				exp = int(math.log10(table[0]['cuts'][cut][0]))
				option = '[table-figures-integer = %d]' % (exp+1)
			file.write('\n\tS%s' % option) # columns for cuts
		file.write('\n\tS[table-figures-integer = 2, table-figures-decimal = 1, table-figures-uncertainty = 2]\n') # column for efficiency
		for FoM in table[0]['results'] :
			file.write('\tS[table-figures-decimal = 1]\n') # columns for results
		file.write('}\n')
		file.write('\t\\toprule\n')

		# first title row
#		file.write('{\\bf Eff.} & ')
		for var in table[0]['cuts'] :
#			file.write('\\multicolumn{2}{c|}{\\bf %s} &' % pl.get_varName(var))
			file.write('\t{\\bf %s} & ' % var)
		file.write('{$\\varepsilon$}')
		for FoM in table[0]['results'] :
			file.write(' & {%s}' % FoM)
		file.write(' \\\\\n')
		file.write('\t\\midrule\n')

		# data rows
		for row in table :
#			file.write('%3d &' % row['ieff'])
			file.write('\t')
			for var in row['cuts'] :
				file.write('%5.0f & ' % row['cuts'][var][0])
			file.write('%4.1f +- %3.1f' % (row['eff'][0]*100., row['eff'][1]*100.))
			for FoM in row['results'] :
				file.write(' & %5.1f' % row['results'][FoM]['6channels'][1])
			file.write(' \\\\\n')
		file.write('\t\\bottomrule\n')
		file.write('\\end{tabular}')


def make_SampleTable(path, samples, name = '') :
	'''list all samples'''

	table_name = '%sSamples.tex' % name
	table_path = path
	if not table_path.endswith('/') : table_path += '/'
	helper.mkdir(table_path)
	pl = ttvStyle.ttvStyle(TeX_switch = True)
	print '[status] writing %s' % table_name
	ngen_digits = int(math.log10(max([s.ngen for s in samples]))) + 1
	xsec_digits = int(math.log10(max([s.xsec for s in samples]))) + 4
	lumi_digits = int(math.log10(max([s.getLumi() for s in samples]))) - 2
	samples.sort(key = lambda s: s.xsec, reverse = True) # sort samples by cross section
	with open(table_path + table_name, 'w') as file :
		timestamp = time.asctime()
		file.write('%!TEX root = ../../Dissertation.tex\n')
#		file.write('\n')
#		file.write(providecommands())
		file.write('\n\n')
		file.write('\\begin{tabular}{\n')
		file.write('\tl\n')
		file.write('\tS[table-number-alignment = center, table-format = %d.2, round-mode = figures, round-precision = 3]\n' % xsec_digits)
		file.write('\tS[table-number-alignment = center, table-format = %d.0]\n' % ngen_digits)
		file.write('\tS[table-number-alignment = center, table-format = %d.0, round-mode = figures, round-precision = 3]\n' % lumi_digits)
		file.write('}\n')
		file.write('\t\\toprule\n')
		file.write('\tSample         & {$\\sigma$ (\\si{\\femto\\barn})} & {$N_{\\text{gen}}$} & {\\intlumi{} (\\si{\\per\\femto\\barn})} \\\\\n')
		file.write('\t\\midrule\n')
		for s in samples :
			file.write('\t%-14s & %28g  &   %15.0f  & %34g  \\\\\n' % (
#				s.name.replace('_', '\\_'),
				pl.get_processCommand(s.name),
				s.xsec * 1000.,
				s.ngen,
				s.getLumi() / 1000.))
		file.write('\t\\bottomrule\n')
		file.write('\\end{tabular}\n')


def providecommands() :
	commands = ''
	commands += '\\providecommand{\\Pe} {\\ensuremath{\\mathrm{e}}}\n'
	commands += '\\providecommand{\\PGm}{\\ensuremath{\\mu}}\n'
	commands += '\\providecommand{\\npp}{\\ensuremath{N_{\\mathrm{pp}}}}\n'
	commands += '\\providecommand{\\npf}{\\ensuremath{N_{\\mathrm{pf}}}}\n'
	commands += '\\providecommand{\\nfp}{\\ensuremath{N_{\\mathrm{fp}}}}\n'
	commands += '\\providecommand{\\nff}{\\ensuremath{N_{\\mathrm{ff}}}}\n'
	commands += '\\providecommand{\\ntt}{\\ensuremath{N_{\\mathrm{tt}}}}\n'
	commands += '\\providecommand{\\ntl}{\\ensuremath{N_{\\mathrm{tl}}}}\n'
	commands += '\\providecommand{\\nlt}{\\ensuremath{N_{\\mathrm{lt}}}}\n'
	commands += '\\providecommand{\\nll}{\\ensuremath{N_{\\mathrm{ll}}}}\n'
	return commands
