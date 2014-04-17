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
		self.do_analysis(effs)
		return

		self.analize_datacards(+1)


	def do_analysis(self, effs) :
		print '[status] setup plotter..'

		# selections
		sels = {}

		# charge mis-ID scale factor
		self.chmid_sf = 1.62
		print 'check: chmid_sf is %f' % self.chmid_sf

		# get fake and prompt ratios
		EWK_SF = {}
		EWK_SF['el']   = self.get_EWK_SF('el')
		EWK_SF['mu17'] = self.get_EWK_SF('mu17')
		EWK_SF['mu24'] = self.get_EWK_SF('mu24')
		self.fpr.fill_ratios(self.get_samples('SingleDoubleMu'), self.get_samples('DoubleEle'), 0, True, EWK_SF)

		for eff in effs :
			cuts = self.read_cuts(self.cutspath + 'cutsGA_Seff%d.txt' % (eff))

			resultspath = self.cutspath + 'results_eff%d.pkl' % eff

			if os.path.exists(resultspath) :
				print '[status] loading results of predictions from %s..' % (resultspath)
				results = helper.load_object(resultspath)

			else :
				print ''
				print '=========================='
				print '| %3.0f%% signal efficiency |' % (eff)
				print '=========================='
				print ''
				print '[status] starting analysis for %2.0f%% signal efficiency point..' % (eff)
				results = {}
				for syst in self.systematics :
					print '[status] making predictions for %s systematic' % (syst)
					self.skim_tree(syst)  # makes sure the skimmed SigEvents tree exists
					systpath = self.path + 'SSDLYields_skim_' + syst + '.root'
					sels[syst] = self.get_selection('eff%d' % eff, cuts, self.systematics[syst])
					results[syst] = self.make_IntPredictions(sels[syst], systpath)
				helper.save_object(results, resultspath)

			# make datacards for each charge-flavor channel
			datacards_6channels = []
			for charge_str in results['Normal'] :
				datacards_3channels = []
				for chan in results['Normal'][charge_str] :
					datacard = self.make_datacard(results, chan, charge_str, 'eff%d' % eff, self.cutspath)
					ch_str   = results['Normal'][charge_str][chan].chan_str
					if charge_str != 'al' and chan != 'al' :
						datacards_6channels.append('%s=%s' % (ch_str, datacard))
					if chan != 'al' :
						datacards_3channels.append('%s=%s' % (ch_str, datacard))
				target_path = '%s/datacard_ssdl_ttW_3channels_%s_eff%s.txt' % (self.cutspath, charge_str, eff)
				self.combine_datacards(datacards_3channels, target_path)
			target_path = '%s/datacard_ssdl_ttW_6channels_eff%s.txt' % (self.cutspath, eff)
			self.combine_datacards(datacards_6channels, target_path)


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
		'''sets selcetion'''

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


	def analize_datacards(self, charge) :

		presel = copy.deepcopy(self.selections['3J1bJ'])
		presel.charge = charge
		presel.sname = 'TTbarW'

		signif_path = self.cutspath + 'signif.pkl'

		if os.path.exists(signif_path) :
			print '[status] loading optimization results..'
			signif = helper.load_object(signif_path)

		else :
			signif = []
			for ieff in range(25, 105, 5) :
				cuts = self.read_cuts(self.cutspath + 'cutsGA_Seff%d.txt' % (ieff))
				sel = self.get_selection('eff%d' % ieff, cuts, 0, True)
				sel.sname = 'TTbarW'
				(eff, eff_err) = self.get_efficiency(self.path + 'SSDLYields_skim_Normal.root', sel, presel)
				print eff
				datacard = self.cutspath + 'datacard_ssdl_ttW_int++_eff%d.txt' % ieff
				signif.append([eff, runCombine.expected_significance(datacard, '')])
			helper.save_object(signif, signif_path)

		self.plot_results(signif)
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
