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
import tables


class optcuts(plotter.plotter) :

	def __init__(self, path, cardfile, cutspath) :
		plotter.plotter.__init__(self, path, cardfile)
		self.cutspath = cutspath
		if not self.cutspath.endswith('/') : self.cutspath += '/'


	def optimize(self, channel) :
		effs = range(25, 105, 5)

		selections = []
		for eff in effs :
			cuts = self.read_cuts(self.cutspath + 'cutsGA_Seff%d.txt' % (eff))
			selections.append((eff, cuts))
		tables.make_CutsTable(self.cutspath, selections)

		self.do_analysis(effs)

		if   channel == 'all'   : charge =  0
		elif channel == 'plus'  : charge = +1
		elif channel == 'minus' : charge = -1

		self.analize_datacards(effs, charge)


	def do_analysis(self, effs) :
		'''makes predictions for all selections from TMVA'''

		print '[status] setup plotter..'

		# charge mis-ID scale factor
		self.chmid_sf = 1.62

		# get fake and prompt ratios
		self.fpr.fill_ratios(self.get_samples('SingleDoubleMu'), self.get_samples('DoubleEle'), 0, True)

		for eff in effs :
			print ''
			print '=========================='
			print '| %3.0f%% signal efficiency |' % (eff)
			print '=========================='
			print ''
			cuts = self.read_cuts(self.cutspath + 'cutsGA_Seff%d.txt' % (eff))
			sel  = self.get_selection('eff%d' % eff, cuts)
			results = self.make_IntPredictions(sel, self.cutspath, '_eff%d' % eff, True)


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
				if sel.charge < 0 : sel.charge = -1  # fixes a bug in cuts files

		return sel


	def analize_datacards(self, effs, charge = 0) :
		'''analyzes int, 3channels and 6channels datacards of ++, -- or without charge selection'''

		presel = copy.deepcopy(self.selections['3J1bJ'])
		presel.charge = charge
		presel.sname = 'TTbarW'

		FoMs = {}
		FoMs['ExpSign'] = {}
		FoMs['XSecErr'] = {}

		for FoM in FoMs :
			signif_path = '%s%s_%s.pkl' % (self.cutspath, FoM, self.get_chargeString(charge))

			if os.path.exists(signif_path) :
				print '[status] loading optimization results..'
				FoMs[FoM] = helper.load_object(signif_path)

			else :
				print ''
				print '========================='
				print '| Analyzing datacards.. |'
				print '========================='
				print ''

				FoMs[FoM]['int'      ] = []
				FoMs[FoM]['3channels'] = []
				if charge == 0 : FoMs[FoM]['6channels'] = []
				for ieff in effs :
					# calculate efficiency of selection
					cuts = self.read_cuts(self.cutspath + 'cutsGA_Seff%d.txt' % (ieff))
					sel = self.get_selection('eff%d' % ieff, cuts, 0, True)
					sel.sname = 'TTbarW'
					(eff, eff_err) = self.get_efficiency(self.path + 'SSDLYields_skim_Normal.root', sel, presel)

					# get expected significance
					for chan in FoMs[FoM] :
						if   chan == '3channels' : charge_suffix = '_' + self.get_chargeString(charge)
						elif chan == 'int'       : charge_suffix =       self.get_chargeString(charge, 1)
						else                     : charge_suffix = ''
						datacard = self.cutspath + 'datacards_eff%d/datacard_ssdl_ttW_%s%s_eff%d.txt' % (ieff, chan, charge_suffix, ieff)
						if FoM == 'ExpSign' : result = runCombine.expected_significance(datacard, '')
						if FoM == 'XSecErr' : result = max(runCombine.signal_strength(datacard, '')[1:])
						FoMs[FoM][chan].append([eff, result])

				helper.save_object(FoMs[FoM], signif_path)

			self.plot_results(FoM, FoMs[FoM], self.get_chargeString(charge))


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


	def plot_results(self, FoM, results, charge_str) :
		pl = ttvplot.ttvplot(self.cutspath, '2L', cms_label = 3)
		canvas = pl.get_canvas()
		canvas.cd()
		scale = 1.
		if FoM == 'XSecErr' : scale = 232.
		h2_axes_gr = ROOT.TH2D("axes_gr", "", 20, 0., 100., 15, 0., 3.*scale)
#		h2_axes_gr = ROOT.TH2D("axes_gr", "", 20, 0., 100., 12, 0., 2.4)
		h2_axes_gr.Draw()
		h2_axes_gr.GetXaxis().SetTitle('#varepsilon_{Signal} [%]')
		h2_axes_gr.GetYaxis().SetTitle('#sigma_{expected}')
		if FoM == 'XSecErr' :  h2_axes_gr.GetYaxis().SetTitle('Exp. Cross Section Error [fb]')
		h2_axes_gr.GetXaxis().SetTitleOffset(1.25)
		h2_axes_gr.GetYaxis().SetTitleOffset(1.25)
		h2_axes_gr.GetXaxis().SetTitleSize(0.046)
		h2_axes_gr.GetYaxis().SetTitleSize(0.046)
		h2_axes_gr.GetXaxis().SetLabelSize(0.04)
		h2_axes_gr.GetYaxis().SetLabelSize(0.04)

		gr_res = {}
		leg_entries = []
		for chan, res in results.iteritems() :

			gr_res[chan] = ROOT.TGraphErrors(0)

			for i, [eff, result] in enumerate(res) :
				gr_res[chan].SetPoint(i, 100.*eff, result*scale)

			gr_res[chan].SetMarkerSize(1.2)
			gr_res[chan].SetMarkerStyle(21)
			if   chan == '6channels' : color = 35
			if   chan == '3channels' : color = 29
			elif chan == 'int'       : color = 40
			gr_res[chan].SetMarkerColor(color)
			gr_res[chan].Draw('P same')
			leg_entries.append([gr_res[chan], chan, 'p'])
		leg = pl.draw_legend(leg_entries)
		pl.draw_cmsLine()
#		raw_input('ok? ')
#		raw_input('ok? ')

		canvas.Print('%sSeff_vs_%s_%s.pdf' % (self.cutspath, FoM, charge_str))
		canvas.Print('%sSeff_vs_%s_%s.png' % (self.cutspath, FoM, charge_str))
#		canvas.Print('%stmp/test.png' % self.path)


if __name__ == '__main__' :
	args = sys.argv

	if ('--help' in args) or ('-h' in args) or ('-d' not in args) or ('-c' not in args) or ('--cuts' not in args) or ('--channel' not in args) :
		print 'usage: ./optcuts.py -d <INPUTDIR> -c <DATACARD> --cuts <CUTSDIR> --channel <CHANNEL>'
		print ''
		print '       -d Directory in which the SSDLYYields.root output from the SDLDumper is.'
		print ''
		print '       -c Datacard with list of samples and corresponding cross sections.'
		print ''
		print '       --cuts Directory in which the output files of TMVA with selection cuts are. Output files are written in this directory too.'
		print ''
		print '       --channel all, plus or minus'
		print ''
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
		if (channel != 'all') and (channel != 'plus') and (channel != 'minus') :
			print '[ERROR] Invalid channel!'
			sys.exit(1)

	opt = optcuts(path, cardfile, cutspath)
	opt.optimize(channel)