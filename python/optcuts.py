#! /usr/bin/python
import os, sys
import plotter
import ratios
import selection
import helper
import runCombine
import ROOT
import ttvplot
import copy


class optcuts(plotter.plotter) :

	def __init__(self, path, cardfile, cutspath) :
		plotter.plotter.__init__(self, path, cardfile)
		self.cutspath = cutspath + '/'


	def optimize(self, channel) :
		effs = range(25, 105, 5)
#		self.do_analysis(effs)
#		return

		self.analize_datacards(effs, +1)


	def do_analysis(self, effs) :
		'''makes predictions for all selections from TMVA'''

		print '[status] setup plotter..'

		# charge mis-ID scale factor
		self.chmid_sf = 1.62

		# get fake and prompt ratios
		EWK_SF = {}
		EWK_SF['el']   = self.get_EWK_SF('el')
		EWK_SF['mu17'] = self.get_EWK_SF('mu17')
		EWK_SF['mu24'] = self.get_EWK_SF('mu24')
		self.fpr.fill_ratios(self.get_samples('SingleDoubleMu'), self.get_samples('DoubleEle'), 0, True, EWK_SF)

		for eff in effs :
			print ''
			print '=========================='
			print '| %3.0f%% signal efficiency |' % (eff)
			print '=========================='
			print ''
			cuts = self.read_cuts(self.cutspath + 'cutsGA_Seff%d.txt' % (eff))
			sel  = self.get_selection('eff%d' % eff, cuts)
			results = self.make_IntPredictions(sel, self.cutspath, '_eff%d' % eff)


	def read_cuts(self, cutfile) :
		'''reads file with cuts'''
		if not os.path.exists(cutfile) :
			print '[ERROR] %s does not exist!' % (cutfile)
			sys.exit(1)
		print '[status] reading cuts from %s' % (cutfile)
		cuts = {}
		with open(cutfile, 'r') as file :
			for i, line in enumerate(file.readlines()) :
				if line[0] is '#' : continue
				splitline = line.split()
				if (len(splitline) != 3) :
					print '[ERROR] check format of %s (line %d)!' % (cutfile, i+1)
					sys.exit(1)
				[var, min, max] = splitline
				min = float(min)
				max = float(max)
				cuts[var] = [min, max]
		return cuts


	def get_selection(self, name, cuts, systflag = 0, charge_sel = False) :

		sel = selection.selection(name)
		sel.systflag = systflag

		for var, cut in cuts.iteritems() :
			if var == 'HT' :
				sel.minHT = float(cut[0])
				sel.maxHT = float(cut[1])
			if var == 'pT1' :
				sel.minPt1 = float(cut[0])
				sel.maxPt1 = float(cut[1])
			if var == 'pT2' :
				sel.minPt2 = float(cut[0])
				sel.maxPt2 = float(cut[1])
			if var == 'NJ' :
				sel.minNjets = int(cut[0])
				sel.maxNjets = int(cut[1])
			if var == 'NbJmed' :
				sel.minNbjetsM = int(cut[0])
				sel.maxNbjetsM = int(cut[1])
			if var == 'Charge' and charge_sel :
				sel.charge = int(cut[0])

		return sel


	def analize_datacards(self, effs, charge = 0) :
		'''analyzes 3channel datacards of ++, -- or without charge selection'''

		presel = copy.deepcopy(self.selections['3J1bJ'])
		presel.charge = charge
		presel.sname = 'TTbarW'

		signif_path = '%ssignif_3channels_%s.pkl' % (self.cutspath, self.get_chargeString(charge))

		if os.path.exists(signif_path) :
			print '[status] loading optimization results..'
			signif = helper.load_object(signif_path)

		else :
			signif = {}
			signif['int'      ] = []
			signif['3channels'] = []
			for ieff in effs :
				# calculate efficiency of selection
				cuts = self.read_cuts(self.cutspath + 'cutsGA_Seff%d.txt' % (ieff))
				sel = self.get_selection('eff%d' % ieff, cuts, 0, True)
				sel.sname = 'TTbarW'
				(eff, eff_err) = self.get_efficiency(self.path + 'SSDLYields_skim_Normal.root', sel, presel)

				# get expected significance
				for chan in signif :
					if   chan == '3channels' : charge_suffix = '_' + self.get_chargeString(charge)
					elif chan == 'int'       : charge_suffix =       self.get_chargeString(charge, 1)
					datacard = self.cutspath + 'datacards_eff%d/datacard_ssdl_ttW_%s%s_eff%d.txt' % (ieff, chan, charge_suffix, ieff)
					signif[chan].append([eff, runCombine.expected_significance(datacard, '')])
			helper.save_object(signif, signif_path)

		self.plot_results(signif['3channels'])
		print '\n\n'
		print signif


	def get_efficiency(self, path, sel, base_sel) :
		'''get efficiency'''

		print '[status] open SigEvents tree from %s' % (path)
		file = ROOT.TFile.Open(path, 'READ')
		tree = file.Get('SigEvents')
		n_after  = tree.GetEntries(sel.get_selectionString())
		n_before = tree.GetEntries(base_sel.get_selectionString())
		ratio = helper.ratioWithBinomErrors(n_after, n_before)
		print '         selection      %10s: %d' % (sel.name, n_after)
		print '         base selection %10s: %d' % (base_sel.name, n_before)
		print '         efficiency: %4.2f +/- %4.2f' % ratio
		return ratio


	def plot_results(self, results) :
		gr_res = ROOT.TGraphErrors(0)
		h2_axes_gr = ROOT.TH2D("axes_gr", "", 20, 0., 100., 25, 0., 5.)
#		h2_axes_gr = ROOT.TH2D("axes_gr", "", 20, 0., 100., 12, 0., 2.4)

		for i, [eff, result] in enumerate(results) :
			print eff, result
			gr_res.SetPoint(i, 100.*eff, result)

		pl = ttvplot.ttvplot(self.cutspath, '2L', cms_label = 3)
		canvas = pl.get_canvas()
#		canvas.cd()

		gr_res.SetMarkerSize(2.)
		gr_res.SetMarkerStyle(21)
		gr_res.SetMarkerColor(29)
		h2_axes_gr.Draw()
		gr_res.Draw('P same')
		pl.draw_cmsLine()
		raw_input('ok? ')
		raw_input('ok? ')

		canvas.Print('%stest.pdf' % self.cutspath)


if __name__ == '__main__' :
	args = sys.argv

	if ('--help' in args) or ('-h' in args) or ('-d' not in args) or ('-c' not in args) or ('--cuts' not in args) or ('--channel' not in args) :
		print 'usage: ./optcuts.py -d <INPUTDIR> -c <DATACARD> -s <CUTSDIR>'
		print ''
		print '       -d Directory in which the SSDLYYields.root output from the SDLDumper is.'
		print ''
		print '       -c Datacard with list of samples and corresponding cross sections.'
		print ''
		print '       -cutsDir Directory in which the output files of TMVA with selection cuts are. Output files are written in this directory too.'
		sys.exit(1)

	if ('-d' in args) and (args[args.index('-d')+1] != '') :
		path = str(args[args.index('-d')+1])
		print path

	if ('-c' in args) and (args[args.index('-c')+1] != '') :
		cardfile = str(args[args.index('-c')+1])
		print cardfile

	if ('--cuts' in args) and (args[args.index('--cuts')+1] != '') :
		cutspath = str(args[args.index('--cuts')+1])

	if ('--channel' in args) and (args[args.index('--channel')+1] != '') :
		channel = str(args[args.index('--channel')+1])

	opt = optcuts(path, cardfile, cutspath)
	opt.optimize(channel)
