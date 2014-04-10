#! /usr/bin/python
#import os, sys, commands, subprocess, math
import ROOT
import sys
import math
import selection
import sample
import helper
import ratios
import result
import time
import copytree
import os
import numpy as np
import ttvplot
import config
import tables


class plotter :
	'''the plotter reads sigtree and produces plots'''

	def __init__(self, path) :
		print '[status] initialize plotter..'
		self.path = path + '/'
		self.ssdlfile = ROOT.TFile.Open(path + '/SSDLYields.root', 'READ')
		if self.ssdlfile == None :
			sys.exit(1)
		self.sigtree = self.ssdlfile.Get('SigEvents')
		print '[status] loaded SigEventsTree with %d events' % (self.sigtree.GetEntries())

		ROOT.gSystem.Load('./FakeRatios.so')

		# samples
		self.samples = {}

		# selections
		self.selections = {}

		# ratios
#		self.fpr = ratios.ratios(path)

		# lumi norm
		self.lumi = 19466.
		self.lumi_HLTMu17       = 24.9
		self.lumi_HLTMu24Eta2p1 = 116.
		self.lumi_HLTEl17Jet30  = 23.845

		# systematic uncertainties
		# TODO: read these from file
		self.RareESyst   = 0.5
		self.RareESyst2  = self.RareESyst * self.RareESyst
		self.FakeESyst   = 0.5
		self.FakeESyst2  = self.FakeESyst * self.FakeESyst
		self.WZESyst     = 0.15
		self.WZESyst2    = self.WZESyst * self.WZESyst
		self.TTZESyst    = 0.5
		self.TTZESyst2   = self.TTZESyst * self.TTZESyst
		self.TTWESyst    = 0.5
		self.TTWESyst2   = self.TTWESyst * self.TTWESyst
		self.ChMisESyst  = 0.3
		self.ChMisESyst2 = self.ChMisESyst * self.ChMisESyst

		# systematics
		self.systematics = {}
		self.systematics['Normal']       = 0
		self.systematics['JetUp']        = 1
		self.systematics['JetDown']      = 2
		self.systematics['JetSmear']     = 3
		self.systematics['BUp']          = 4
		self.systematics['BDown']        = 5
		self.systematics['LepUp']        = 6
		self.systematics['LepDown']      = 7
#		self.systematics['METUp']        = 8
#		self.systematics['METDown']      = 9
		self.systematics['PileupUp']     = 10
		self.systematics['PileupDown']   = 11
		self.systematics['JetSmearUp']   = 12
		self.systematics['JetSmearDown'] = 13
		self.systematics['MuUp']         = 14
		self.systematics['MuDown']       = 15
		self.systematics['ElUp']         = 16
		self.systematics['ElDown']       = 17

		# systematic uncertainties for datacard
		self.lumi_syst     = 1.026
		self.scale_syst_up = 1.027
		self.scale_syst_dn = 0.980
		self.tmass_syst_up = 1.019
		self.tmass_syst_dn = 0.983
		self.ltrig_syst    = 1.030
		self.pdf_syst      = 1.015
		self.gen_syst      = 1.050
		self.wz_pdf_syst   = 1.041

		# charge strings
		self.charges = {}
		for charge in range(-1, 2) :
			self.charges[self.get_chargeString(charge)] = charge


	def do_analysis(self, cardfile) :
		print '[status] starting analysis..'

		# samples
		self.samples = self.readDatacard(cardfile)
		self.read_ngen()

		# selections
		sels = {}
#		sels['1J0bJ'    ] = selection.selection(name = '1J0bJ'  , minNjets = 1)
#		sels['2J0bJ'    ] = selection.selection(name = '2J0bJ'  , minNjets = 2)  # loose selection
#		sels['2JnobJ_ee'] = selection.selection(name = '2JnobJ' , minNjets = 2, maxNbjetsM = 0, flavor = 2)
		sels['3J1bJ'    ] = selection.selection(name = '3J1bJ'  , minNjets = 3, minNbjetsM = 1)  # pre-selection
#		sels['3J1bJ_ee' ] = selection.selection(name = '3J1bJ'  , minNjets = 3, minNbjetsM = 1, flavor = 2)
		sels['final'    ] = selection.selection(name = 'final'  , minNjets = 3, minNbjetsL = 1, minNbjetsM = 1, minPt1 = 40., minPt2 = 40., minHT = 155., charge = 0)
#		sels['final++'  ] = selection.selection(name = 'final++', minNjets = 3, minNbjetsL = 1, minNbjetsM = 1, minPt1 = 40., minPt2 = 40., minHT = 155., charge = +1)
#		sels['final--'  ] = selection.selection(name = 'final--', minNjets = 3, minNbjetsL = 1, minNbjetsM = 1, minPt1 = 40., minPt2 = 40., minHT = 155., charge = -1)
		#print presel.get_selectionString()
#		self.selections[presel.name] = presel

		self.chmid_sf = 1.62
#		self.chmid_sf = 1.

		EWK_SF = {}
		EWK_SF['el']   = self.get_EWK_SF('el')
		EWK_SF['mu17'] = self.get_EWK_SF('mu17')
		EWK_SF['mu24'] = self.get_EWK_SF('mu24')

		applyEwkSubtr = True
#		self.fill_ratios(self.get_samples('SingleDoubleMu'), self.get_samples('DoubleEle'), 0, True, EWK_SF)
		self.fpr = ratios.ratios(self.path, self.samples)
		self.fpr.fill_ratios(self.get_samples('SingleDoubleMu'), self.get_samples('DoubleEle'), 0, True, EWK_SF)

		print 'MufRatio: %f +/- %f' % (self.fpr.MufRatio, self.fpr.MufRatioE)
		print 'ElfRatio: %f +/- %f' % (self.fpr.ElfRatio, self.fpr.ElfRatioE)
		print 'MupRatio: %f +/- %f' % (self.fpr.MupRatio, self.fpr.MupRatioE)
		print 'ElpRatio: %f +/- %f' % (self.fpr.ElpRatio, self.fpr.ElpRatioE)

#		c1 = ROOT.TCanvas("canvas", "canvas", 0, 0, 800, 800)
#		c1.Divide(2, 2)
#		c1.cd(1)
#		self.h2_MufRatio.Draw('colztext')
#		c1.cd(2)
#		self.h2_ElfRatio.Draw('colztext')
#		c1.cd(3)
#		self.h2_MupRatio.Draw('colztext')
#		c1.cd(4)
#		self.h2_ElpRatio.Draw('colztext')
#


#		# produce results tree
		self.skim_tree('Normal')  # makes sure the skimmed SigEvents tree exists
		restree_path = {}
		restree_path['1J0bJ'] = self.path + 'SSDLResults.root'
		restree_path['2J0bJ'] = self.path + 'SSDLResults_2J0bJ.root'
		restree_path['3J1bJ'] = self.path + 'SSDLResults_3J1bJ.root'
		if not os.path.exists(restree_path['1J0bJ']) :
			self.make_IntPredictions(sels['1J0bJ'], self.path + 'SSDLYields_skim_Normal.root', True)
		if not os.path.exists(restree_path['2J0bJ']) :
			copytree.copytree(restree_path['1J0bJ'], restree_path['2J0bJ'], 'Results', 'NJ > 1')
		if not os.path.exists(restree_path['3J1bJ']) :
			copytree.copytree(restree_path['1J0bJ'], restree_path['3J1bJ'], 'Results', 'NJ > 2 && NbJmed > 0')

		for name, sel in sels.iteritems() :
			res_path = restree_path['3J1bJ']
			if sel.name.startswith('2J') : res_path = restree_path['2J0bJ']
			if sel.name.startswith('1J') : res_path = restree_path['1J0bJ']
			self.plot_predictions(res_path, sel)

		return


#		nosyst_path = self.path + 'SSDLYields_Normal.root'
#		if not os.path.exists(nosyst_path) :
#			copytree.copytree(self.path + 'SSDLYields.root', nosyst_path, 'SigEvents', 'SystFlag == 0 && (SType < 3 || TLCat == 0) && (SType < 3 || Flavor < 3) && (SType < 3 || SType == 15)')
#
#		nosystfile = ROOT.TFile.Open(nosyst_path, 'READ')
#		nosysttree = nosystfile.Get('SigEvents')
#
#		self.make_IntPredictions(loosesel, nosysttree, True)

		#return

#		# plots results from results tree
#		resultfile = ROOT.TFile.Open(self.path + 'SSDLResults.root', 'READ')
#		resulttree = resultfile.Get('Results')
##		self.plot_predictions(resulttree, finalsel)
##		self.plot_predictions(resulttree, presel)
##		self.plot_predictions(resulttree, loosesel)
##		self.plot_predictions(resulttree, presel_ee)
##		self.plot_predictions(resulttree, sels['2JnobJ'])
#		self.plot_predictions(resulttree, sels['final++'])
#		self.plot_predictions(resulttree, sels['final--'])
#		return




		results = {}
		resultspath = self.path + 'IntPredictions/results.pkl'

		sels = {}

		if os.path.exists(resultspath) :
			print '[status] loading results of predictions from %s..' % (resultspath)
			results = helper.load_object(resultspath)

		else :
			for syst in self.systematics :
#				if syst != 'Normal' : continue
				print '[status] making predictions for %s systematic' % (syst)
				self.skim_tree(syst)  # makes sure the skimmed SigEvents tree exists
				systpath = self.path + 'SSDLYields_skim_' + syst + '.root'
				sels[syst] = selection.selection(name = 'final_' + syst, minNjets = 3, minNbjetsL = 1, minNbjetsM = 1, minPt1 = 40., minPt2 = 40., minHT = 155., systflag = self.systematics[syst])
				results[syst] = self.make_IntPredictions(sels[syst], systpath)

			helper.save_object(results, resultspath)

		# make datacards for each charge-flavor channel
		for ch_str in results['Normal'] :
			for chan in results['Normal'][ch_str] :
				self.make_datacard(results, chan, ch_str)

		# make table of observation and predictions
		tables.make_ObsPredTable(self.path, results['Normal'])

		# make overview table of systematic studies
		tables.make_SystTable(self.path, results, 'al', 'al')

		print '++:   Observed: %3d   Predicted: %5.1f +/- %4.1f' % (results['Normal']['++']['al'].obs, results['Normal']['++']['al'].tot + results['Normal']['++']['al'].ttw, math.sqrt(results['Normal']['++']['al'].tot_err*results['Normal']['++']['al'].tot_err + results['Normal']['++']['al'].ttw_err*results['Normal']['++']['al'].ttw_err))
		print '--:   Observed: %3d   Predicted: %5.1f +/- %4.1f' % (results['Normal']['--']['al'].obs, results['Normal']['--']['al'].tot + results['Normal']['--']['al'].ttw, math.sqrt(results['Normal']['--']['al'].tot_err*results['Normal']['--']['al'].tot_err + results['Normal']['--']['al'].ttw_err*results['Normal']['--']['al'].ttw_err))

#
#		raw_input('ok? ')


	def skim_tree(self, syst = '') :
		'''skim SigEvents tree'''

		skimtree_path = self.path + 'SSDLYields_skim.root'
		if not os.path.exists(skimtree_path) :
			print '[status] creating skimmed tree for di-boson, rare MC (SType 15) and data..'
			print '         Data: same-sign, opposite-sign, tight-tight, tight-loose, loose-tight, loose-loose events (only nominal, no systematics)'
			print '         MC:   only same-sign, tight-tight events for all systematics'
			copytree.copytree(self.path + 'SSDLYields.root', skimtree_path, 'SigEvents', '(SType < 3 || SType == 15) && (SType < 3 || Flavor < 3) && (SType < 3 || TLCat == 0) && (SType > 2 || SystFlag == 0)')

		if syst != '' :
			systtree_path = self.path + 'SSDLYields_skim_' + syst + '.root'

			if not os.path.exists(systtree_path) :
				print '[status] creating skimmed tree file for %s systematic..' % (syst)
				copytree.copytree(skimtree_path, systtree_path, 'SigEvents', 'SystFlag == %d' % (self.systematics[syst]))


	def readDatacard(self, cardfile, verbose = 0) :
		'''reads datacard file'''
		print '[status] reading datacard %s' % (cardfile)
		samples = {}
		## line format:
		##
		## data:
		##  0: sample name
		##  1: input file
		##  2: data mc type
		##  3: channel (0: DoubleMu, 1: DoubleEle, 2: MuEG, 5: SingleMu)
		##
		## mc:
		##  0: sample name
		##  1: input file
		##  2: data mc type
		##  3: channel (-1)
		##  4: cross section
		with open(cardfile, 'r') as file :
			for i, line in enumerate(file.readlines()) :
				if line[0] is '#' : continue
				splitline = line.split()
				if (len(splitline) < 4) or (len(splitline) > 5) or (len(splitline) == 4 and int(splitline[2]) != 0) or (len(splitline) == 5 and int(splitline[2]) == 0) :
					print '[error] check format of your datacard (line %d)!' % (i)
					sys.exit(1)
				[name, inputfile, datamc, channel] = splitline[:4]
				datamc = int(datamc)
				channel = int(channel)
				xsec = -1.
				if len(splitline) == 5 :
					xsec = float(splitline[4])
				samples[name] = sample.sample(name = name, datamc = datamc, channel = channel, xsec = xsec, ngen = -1)
				if verbose > 0 : print samples[name]
		return samples


	def read_ngen(self) :
		for name, s in self.samples.iteritems() :
			if s.datamc is 0 : continue
			s.ngen = int(self.ssdlfile.Get(s.name+'/'+s.name+'_EventCount').GetEntries())


	def get_samples(self, channel) :
		samplelist = []
		if channel == 'DoubleMu' :
			for name, sample in self.samples.iteritems() :
				if (sample.datamc == 0) and (sample.channel == 0) : samplelist.append(sample.name)
		if channel == 'DoubleEle' :
			for name, sample in self.samples.iteritems() :
				if (sample.datamc == 0) and (sample.channel == 1) : samplelist.append(sample.name)
		if channel == 'MuEG' :
			for name, sample in self.samples.iteritems() :
				if (sample.datamc == 0) and (sample.channel == 2) : samplelist.append(sample.name)
		if channel == 'SingleMu' :
			for name, sample in self.samples.iteritems() :
				if (sample.datamc == 0) and (sample.channel == 5) : samplelist.append(sample.name)
		if channel == 'SingleDoubleMu' :
			for name, sample in self.samples.iteritems() :
				if (sample.datamc == 0) and ((sample.channel == 0) or (sample.channel == 5)) : samplelist.append(sample.name)
		return samplelist


	def get_channelString(self, flavor, opt = 0) :
		if opt is 0 :
			if flavor is 0 : return 'Muon'
			if flavor is 1 : return 'ElMu'
			if flavor is 2 : return 'Elec'
			return ''

		if opt is 1 :
			if flavor is 0 : return 'MM'
			if flavor is 1 : return 'EM'
			if flavor is 2 : return 'EE'
			if flavor is 3 : return 'MM_OS'
			if flavor is 4 : return 'EM_OS'
			if flavor is 5 : return 'EE_OS'
			return ''


	def get_chargeString(self, charge, opt = 0) :
		if opt is 0:
			if charge ==  0 : return 'al'
			if charge == +1 : return '++'
			if charge == -1 : return '--'
			return ''


	def read_histos(self) :
		'''reads histograms for all samples from SSDLYields.root'''
		print '[status] reading histograms..'


	def book_resultsTree(self) :
		# create branches
		self.ResTree_Weight    = np.zeros(1, dtype=float); self.results_tree.Branch('Weight',      self.ResTree_Weight     , 'Weight/D'   );
		self.ResTree_ObsPred   = np.zeros(1, dtype=int  ); self.results_tree.Branch('ObsPred',     self.ResTree_ObsPred    , 'ObsPred/I'  );
		self.ResTree_SystFlag  = np.zeros(1, dtype=int  ); self.results_tree.Branch('SystFlag',    self.ResTree_SystFlag   , 'SystFlag/I' );
		self.ResTree_SName     = ROOT.std.string(       ); self.results_tree.Branch('SName',       self.ResTree_SName);
		self.ResTree_SType     = np.zeros(1, dtype=int  ); self.results_tree.Branch('SType',       self.ResTree_SType      , 'SType/I'    );
		self.ResTree_Run       = np.zeros(1, dtype=int  ); self.results_tree.Branch('Run',         self.ResTree_Run        , 'Run/I'      );
		self.ResTree_LS        = np.zeros(1, dtype=int  ); self.results_tree.Branch('LS',          self.ResTree_LS         , 'LS/I'       );
		self.ResTree_Event     = np.zeros(1, dtype=int  ); self.results_tree.Branch('Event',       self.ResTree_Event      , 'Event/I'    );
		self.ResTree_Flavor    = np.zeros(1, dtype=int  ); self.results_tree.Branch('Flavor',      self.ResTree_Flavor     , 'Flavor/I'   );
		self.ResTree_Charge    = np.zeros(1, dtype=int  ); self.results_tree.Branch('Charge',      self.ResTree_Charge     , 'Charge/I'   );
		self.ResTree_CFChan    = np.zeros(1, dtype=int  ); self.results_tree.Branch('CFChan',      self.ResTree_CFChan     , 'CFChan/I'   );
		self.ResTree_TLCat     = np.zeros(1, dtype=int  ); self.results_tree.Branch('TLCat',       self.ResTree_TLCat      , 'TLCat/I'    );
#		self.ResTree_ZVeto     = np.zeros(1, dtype=int  ); self.results_tree.Branch('PassZVeto',   self.ResTree_ZVeto      , 'PassZVeto/I');
		self.ResTree_HT        = np.zeros(1, dtype=float); self.results_tree.Branch('HT',          self.ResTree_HT         , 'HT/D'       );
		self.ResTree_MET       = np.zeros(1, dtype=float); self.results_tree.Branch('MET',         self.ResTree_MET        , 'MET/D'      );
		self.ResTree_NJ        = np.zeros(1, dtype=int  ); self.results_tree.Branch('NJ',          self.ResTree_NJ         , 'NJ/I'       );
		self.ResTree_NbJ       = np.zeros(1, dtype=int  ); self.results_tree.Branch('NbJ',         self.ResTree_NbJ        , 'NbJ/I'      );
		self.ResTree_NbJmed    = np.zeros(1, dtype=int  ); self.results_tree.Branch('NbJmed',      self.ResTree_NbJmed     , 'NbJmed/I'   );
		self.ResTree_Mll       = np.zeros(1, dtype=float); self.results_tree.Branch('Mll',         self.ResTree_Mll        , 'Mll/D'      );
		self.ResTree_pT1       = np.zeros(1, dtype=float); self.results_tree.Branch('pT1',         self.ResTree_pT1        , 'pT1/D'      );
		self.ResTree_pT2       = np.zeros(1, dtype=float); self.results_tree.Branch('pT2',         self.ResTree_pT2        , 'pT2/D'      );
#		self.ResTree_PFIso1    = np.zeros(1, dtype=float); self.results_tree.Branch('PFIso1',      self.ResTree_PFIso1     , 'PFIso1/D'   );
#		self.ResTree_PFIso2    = np.zeros(1, dtype=float); self.results_tree.Branch('PFIso2',      self.ResTree_PFIso2     , 'PFIso2/D'   );
#		self.ResTree_D01       = np.zeros(1, dtype=float); self.results_tree.Branch('D01',         self.ResTree_D01        , 'D01/D'      );
#		self.ResTree_D02       = np.zeros(1, dtype=float); self.results_tree.Branch('D02',         self.ResTree_D02        , 'D02/D'      );
#		self.ResTree_Rho       = np.zeros(1, dtype=float); self.results_tree.Branch('Rho',         self.ResTree_Rho        , 'Rho/D'      );
		self.ResTree_MTLep1    = np.zeros(1, dtype=float); self.results_tree.Branch('MTLep1',      self.ResTree_MTLep1     , 'MTLep1/D'   );
		self.ResTree_MTLep2    = np.zeros(1, dtype=float); self.results_tree.Branch('MTLep2',      self.ResTree_MTLep2     , 'MTLep2/D'   );
#		self.ResTree_dRbJl1    = np.zeros(1, dtype=float); self.results_tree.Branch('dRbJl1',      self.ResTree_dRbJl1     , 'dRbJl1/D'   );
#		self.ResTree_dRbJl2    = np.zeros(1, dtype=float); self.results_tree.Branch('dRbJl2',      self.ResTree_dRbJl2     , 'dRbJl2/D'   );
#		self.ResTree_BetaStar1 = np.zeros(1, dtype=float); self.results_tree.Branch('BetaStar1',   self.ResTree_BetaStar1  , 'BetaStar1/D');
#		self.ResTree_BetaStar2 = np.zeros(1, dtype=float); self.results_tree.Branch('BetaStar2',   self.ResTree_BetaStar2  , 'BetaStar2/D');
#		self.ResTree_BetaStar3 = np.zeros(1, dtype=float); self.results_tree.Branch('BetaStar3',   self.ResTree_BetaStar3  , 'BetaStar3/D');
#		self.ResTree_BetaStar4 = np.zeros(1, dtype=float); self.results_tree.Branch('BetaStar4',   self.ResTree_BetaStar4  , 'BetaStar4/D');
#		self.ResTree_BetaStar5 = np.zeros(1, dtype=float); self.results_tree.Branch('BetaStar5',   self.ResTree_BetaStar5  , 'BetaStar5/D');
		self.ResTree_NVrtx     = np.zeros(1, dtype=int  ); self.results_tree.Branch('NVrtx',       self.ResTree_NVrtx      , 'NVrtx/I'    );
		self.ResTree_M3        = np.zeros(1, dtype=float); self.results_tree.Branch('M3',          self.ResTree_M3         , 'M3/D'       );


	def get_EWK_SF(self, chan_str) :
		print '[status] calculating EWK scale factor for %s trigger' % (chan_str)
		samples_wjets = []
		samples_zjets = []
		samples_qcd = [] # TODO
		if chan_str is 'el'   : samples_data = self.get_samples('DoubleEle')
		if chan_str is 'mu17' : samples_data = self.get_samples('DoubleMu')
		if chan_str is 'mu24' : samples_data = self.get_samples('SingleMu')
		samples_wjets.append('WJets')
		samples_zjets.append('DYJets')

		(h_ntight_data , h_nloose_data ) = self.get_fRatioPlots(samples_data , chan_str, 'MT_MET30')
		(h_ntight_wjets, h_nloose_wjets) = self.get_fRatioPlots(samples_wjets, chan_str, 'MT_MET30')
		(h_ntight_zjets, h_nloose_zjets) = self.get_fRatioPlots(samples_zjets, chan_str, 'MT_MET30')

		bin_min = h_ntight_data.FindBin(60.)
		bin_max = h_ntight_data.FindBin(90.)-1

		n_data = h_ntight_data.Integral(bin_min, bin_max)
		n_mc   = h_ntight_wjets.Integral(bin_min, bin_max) + h_ntight_zjets.Integral(bin_min, bin_max)

		print '         %s SF = data / (WJets + DYJets) = %f / %f = %f' % (chan_str, n_data, n_mc, n_data/n_mc)

		return n_data / n_mc


	def get_fRatioPlots(self, samples, chan_str, ratiovar) :
		'''gets ntight and loose histograms for various variables'''
		## chan_str = 'mu', 'el', 'mu17', 'mu24'
		lumi = self.lumi
		if chan_str is 'el'   : lumi = self.lumi_HLTEl17Jet30
		if chan_str is 'mu17' : lumi = self.lumi_HLTMu17
		if chan_str is 'mu24' : lumi = self.lumi_HLTMu24Eta2p1
		for i, s in enumerate(samples) :
			scale = lumi / self.samples[s].getLumi()
			if self.samples[s].datamc == 0 : scale = 1.
			if i is 0 :
				h_ntight = self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_ntight_'+ratiovar).Clone()
				h_nloose = self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_nloose_'+ratiovar).Clone()
				h_ntight.Scale(scale)
				h_nloose.Scale(scale)
			else :
				h_ntight.Add(self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_ntight_'+ratiovar), scale)
				h_nloose.Add(self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_nloose_'+ratiovar), scale)
		return (h_ntight, h_nloose)


	def get_ChMisID_SF(self) :
		foo = 0


	def calculateChMisIdProb(self, samples, chmid_reg, chmid_sf = 1.) :
		ospair = 0.
		sspair = 0.
		ospair_e = 0.
		sspair_e = 0.

		for s in samples :
			if self.samples[s].datamc is 0 : scale = 1.
			else                           : scale = self.lumi / self.samples[s].getLumi()

			ospairstmp = self.ssdlfile.Get(s+'/ChMisID/'+s+'_EE_ospairs').Clone()
			sspairstmp = self.ssdlfile.Get(s+'/ChMisID/'+s+'_EE_sspairs').Clone()

			if chmid_reg is 'BB' :
				ospair   += ospairstmp.GetBinContent(1,1) * scale
				sspair   += sspairstmp.GetBinContent(1,1) * scale
				ospair_e += ospairstmp.GetBinError(1,1)   * scale
				sspair_e += sspairstmp.GetBinError(1,1)   * scale

			elif chmid_reg is 'EB' :
				ospair   += (ospairstmp.GetBinContent(1,2) + ospairstmp.GetBinContent(2,1)) * scale
				sspair   += (sspairstmp.GetBinContent(1,2) + sspairstmp.GetBinContent(2,1)) * scale
				ospair_e += (ospairstmp.GetBinError(1,2)   + ospairstmp.GetBinError(2,1)  ) * scale
				sspair_e += (sspairstmp.GetBinError(1,2)   + sspairstmp.GetBinError(2,1)  ) * scale

			elif chmid_reg is 'EE' :
				ospair   += ospairstmp.GetBinContent(2,2) * scale
				sspair   += sspairstmp.GetBinContent(2,2) * scale
				ospair_e += ospairstmp.GetBinError(2,2)   * scale
				sspair_e += sspairstmp.GetBinError(2,2)   * scale

		(chmid, chmide) = helper.ratioWithBinomErrors(sspair, ospair)

		# Divide to get the per-electron probability...
		chmid  = chmid  / 2.
		chmide = chmide / 2.

		# apply correction scale factor
		chmid  = chmid  * chmid_sf
		chmide = chmide * chmid_sf

		return (chmid, chmide)


	def makeRatioControlPlots(self) :
		foo = 0


	def make_IntPredictions(self, sel, treepath, write_ResTree = False) :
		'''oberservation and prediction for different selections'''

		print '[status] getting observations and predictions for results tree with selection %s:' % (sel.name)
		print sel

		print '[status] open SigEvents tree from %s' % (treepath)
		sigfile = ROOT.TFile.Open(treepath, 'READ')
		sigtree = sigfile.Get('SigEvents')
		print '[status] loaded %s tree with %d events' % (sigtree.GetName(), sigtree.GetEntries())

		##################
		# INIT VARIABLES #
		##################

		FR = ROOT.FakeRatios()

		res = {}
		for ch_str, charge in self.charges.iteritems() :
			res[ch_str] = {}
			charge_str = ''
			if charge > 0 : charge_str = '+'
			if charge < 0 : charge_str = '-'
			res[ch_str]['al'] = result.result('al', charge, 'int'+charge_str+charge_str)
			res[ch_str]['mm'] = result.result('mm', charge, 'm'+charge_str+'m'+charge_str)
			res[ch_str]['em'] = result.result('em', charge, 'e'+charge_str+'m'+charge_str)
			res[ch_str]['ee'] = result.result('ee', charge, 'e'+charge_str+'e'+charge_str)

		# all rares, ttW and ttZ yields (tight-tight)
		rares = {}
		rares_npass = {}

		#######################
		# ChMisID predictions #
		#######################

		(fbb, fbbE) = self.calculateChMisIdProb(self.get_samples('DoubleEle'), 'BB', self.chmid_sf)
		(feb, febE) = self.calculateChMisIdProb(self.get_samples('DoubleEle'), 'EB', self.chmid_sf)
		(fee, feeE) = self.calculateChMisIdProb(self.get_samples('DoubleEle'), 'EE', self.chmid_sf)

#		(fbb_mc, fbbE_mc) = calculateChMisIdProb(fMCBG, 'BB')
#		(feb_mc, febE_mc) = calculateChMisIdProb(fMCBG, 'EB')
#		(fee_mc, feeE_mc) = calculateChMisIdProb(fMCBG, 'EE')

		if sel.systflag is 0 and not write_ResTree :
			print '[status] getting opposite-sign yields..'

			# EM OS
			nt2_em_BB_os = sigtree.GetEntries(sel.get_selectionString((4,0)))
			nt2_em_EE_os = sigtree.GetEntries(sel.get_selectionString((4,1)))

			# EE OS
			nt2_ee_BB_os = sigtree.GetEntries(sel.get_selectionString((5,0)))
			nt2_ee_EB_os = sigtree.GetEntries(sel.get_selectionString((5,4)))
			nt2_ee_EE_os = sigtree.GetEntries(sel.get_selectionString((5,3)))

			for ch_str, charge in self.charges.iteritems() :

				# charge factor: only takes half of the ChMisID prediction if a charge selection is applied
				chargeFactor = 1.
				if charge != 0: chargeFactor = 0.5

				# Simple error propagation assuming error on number of events is sqrt(N)
				nt2_ee_chmid    = chargeFactor * (2*fbb* nt2_ee_BB_os                           + 2*fee*nt2_ee_EE_os                      + 2*feb*nt2_ee_EB_os)
				nt2_ee_chmid_e1 = chargeFactor * math.sqrt( 4*fbb*fbb * FR.getEStat2(nt2_ee_BB_os)  + 4*fee*fee * FR.getEStat2(nt2_ee_EE_os)    + 4*feb*feb * FR.getEStat2(nt2_ee_EB_os) )  # stat only
				nt2_ee_chmid_e2 = chargeFactor * math.sqrt( 4*fbbE*fbbE * nt2_ee_BB_os*nt2_ee_BB_os + 4*feeE*feeE * nt2_ee_EE_os*nt2_ee_EE_os   + 4*febE*febE * nt2_ee_EB_os*nt2_ee_EB_os + self.ChMisESyst2 * nt2_ee_chmid*nt2_ee_chmid/(chargeFactor*chargeFactor) )  # syst only

				nt2_em_chmid    = chargeFactor * (fbb * nt2_em_BB_os + fee * nt2_em_EE_os)
				nt2_em_chmid_e1 = chargeFactor * math.sqrt( fbb*fbb * FR.getEStat2(nt2_em_BB_os)  + fee*fee*FR.getEStat2(nt2_em_EE_os) )
				nt2_em_chmid_e2 = chargeFactor * math.sqrt( fbbE*fbbE * nt2_em_BB_os*nt2_em_BB_os + feeE*feeE * nt2_em_EE_os*nt2_em_EE_os + 0.25*self.ChMisESyst2*nt2_em_chmid*nt2_em_chmid/(chargeFactor*chargeFactor) )

				res[ch_str]['al'].cmid     = nt2_ee_chmid + nt2_em_chmid
				res[ch_str]['em'].cmid     = nt2_em_chmid
				res[ch_str]['ee'].cmid     = nt2_ee_chmid
				res[ch_str]['al'].cmid_err = math.sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2)
				res[ch_str]['em'].cmid_err = math.sqrt(nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2)
				res[ch_str]['ee'].cmid_err = math.sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2)

		else :
			print '[status] skipped charge mis-ID prediction for systematics'

		################
		# Results tree #
		################

		if write_ResTree :
			self.results_file = ROOT.TFile(self.path + 'SSDLResults.root', 'RECREATE')
			self.results_tree = ROOT.TTree('Results', 'ResultsTree')
			self.book_resultsTree()

		####################################
		# Loop over skimmed SigEvents tree #
		####################################

		last_sample = ''
		nevents = {}
		for s in self.samples : nevents[s] = 0

		# remove this again
		HLTSF = 0.

		print '[status] looping over skimmed tree..'
#		for i, event in enumerate(self.skimtree) :
		for i, event in enumerate(sigtree) :
			HLTSF = event.HLTSF
			nevents[str(event.SName)] += 1
			if last_sample != str(event.SName) :
				print '[status] processing %s..' % (event.SName)
				last_sample = str(event.SName)

			if not sel.passes_selection(event = event, ttLeptons = False, OSwoZVeto = True) : continue
			chan = self.get_channelString(int(event.Flavor))

			# SET VARIABLES FOR RESULTS TREE
			if write_ResTree :
				self.ResTree_SystFlag [0] = event.SystFlag
				self.ResTree_SName.assign(event.SName)
				self.ResTree_SType    [0] = event.SType
				self.ResTree_Run      [0] = event.Run
				self.ResTree_LS       [0] = event.LS
				self.ResTree_Event    [0] = event.Event
				if event.Flavor < 3 :
					self.ResTree_Flavor   [0] = event.Flavor
				else :
					self.ResTree_Flavor   [0] = event.Flavor - 3
				self.ResTree_Charge   [0] = event.Charge
				if event.Charge > 0 :
					self.ResTree_CFChan   [0] = self.ResTree_Flavor[0]
				else :
					self.ResTree_CFChan   [0] = self.ResTree_Flavor[0] + 3
				self.ResTree_TLCat    [0] = event.TLCat
#				self.ResTree_ZVeto    [0] = event.PassZVeto
				self.ResTree_HT       [0] = event.HT
				self.ResTree_MET      [0] = event.MET
				self.ResTree_NJ       [0] = event.NJ
				self.ResTree_NbJ      [0] = event.NbJ
				self.ResTree_NbJmed   [0] = event.NbJmed
				self.ResTree_Mll      [0] = event.Mll
				self.ResTree_pT1      [0] = event.pT1
				self.ResTree_pT2      [0] = event.pT2
#				self.ResTree_PFIso1   [0] = event.PFIso1
#				self.ResTree_PFIso2   [0] = event.PFIso2
#				self.ResTree_D01      [0] = event.D01
#				self.ResTree_D02      [0] = event.D02
#				self.ResTree_Rho      [0] = event.Rho
				self.ResTree_MTLep1   [0] = event.MTLep1
				self.ResTree_MTLep2   [0] = event.MTLep2
#				self.ResTree_dRbJl1   [0] = event.dRbJl1
#				self.ResTree_dRbJl2   [0] = event.dRbJl2
#				self.ResTree_BetaStar1[0] = event.BetaStar1
#				self.ResTree_BetaStar2[0] = event.BetaStar2
#				self.ResTree_BetaStar3[0] = event.BetaStar3
#				self.ResTree_BetaStar4[0] = event.BetaStar4
#				self.ResTree_BetaStar5[0] = event.BetaStar5
				self.ResTree_NVrtx    [0] = event.NVrtx
				self.ResTree_M3       [0] = event.M3v2

			# GET ALL DATA EVENTS
			if event.SType < 3 :
				if event.Flavor < 3 :

					# fake, prompt predictions
					if chan is 'ElMu' :
						f1 = self.fpr.get_fRatio('Muon', event.pT1, event.eta1, 0)
						f2 = self.fpr.get_fRatio('Elec', event.pT2, event.eta2, 0)
						p1 = self.fpr.get_pRatio('Muon', event.pT1, 0)
						p2 = self.fpr.get_pRatio('Elec', event.pT2, 0)
					else :
						f1 = self.fpr.get_fRatio(chan, event.pT1, event.eta1, 0)
						f2 = self.fpr.get_fRatio(chan, event.pT2, event.eta2, 0)
						p1 = self.fpr.get_pRatio(chan, event.pT1, 0)
						p2 = self.fpr.get_pRatio(chan, event.pT2, 0)

					npp = FR.getWpp(event.TLCat, f1, f2, p1, p2)
					npf = FR.getWpf(event.TLCat, f1, f2, p1, p2)
					nfp = FR.getWfp(event.TLCat, f1, f2, p1, p2)
					nff = FR.getWff(event.TLCat, f1, f2, p1, p2)

					for ch_str, charge in self.charges.iteritems() :
						if charge != 0 and event.Charge != charge : continue

						# int
						if event.TLCat is 0 : res[ch_str]['al'].nt2  += 1
						if event.TLCat is 1 : res[ch_str]['al'].nt10 += 1
						if event.TLCat is 2 : res[ch_str]['al'].nt01 += 1
						if event.TLCat is 3 : res[ch_str]['al'].nt0  += 1

						# MM
						if event.Flavor is 0 :
							res[ch_str]['mm'].npp += npp;
							res[ch_str]['mm'].npf += npf;
							res[ch_str]['mm'].nfp += nfp;
							res[ch_str]['mm'].nff += nff;
							if event.TLCat is 0 : res[ch_str]['mm'].nt2  += 1
							if event.TLCat is 1 : res[ch_str]['mm'].nt10 += 1
							if event.TLCat is 2 : res[ch_str]['mm'].nt01 += 1
							if event.TLCat is 3 : res[ch_str]['mm'].nt0  += 1

						# EM
						if event.Flavor is 1 :
							res[ch_str]['em'].npp += npp
							res[ch_str]['em'].npf += npf
							res[ch_str]['em'].nfp += nfp
							res[ch_str]['em'].nff += nff
							if event.TLCat is 0 : res[ch_str]['em'].nt2  += 1
							if event.TLCat is 1 : res[ch_str]['em'].nt10 += 1
							if event.TLCat is 2 : res[ch_str]['em'].nt01 += 1
							if event.TLCat is 3 : res[ch_str]['em'].nt0  += 1

						# EE
						if event.Flavor is 2 :
							res[ch_str]['ee'].npp += npp
							res[ch_str]['ee'].npf += npf
							res[ch_str]['ee'].nfp += nfp
							res[ch_str]['ee'].nff += nff
							if event.TLCat is 0 : res[ch_str]['ee'].nt2  += 1
							if event.TLCat is 1 : res[ch_str]['ee'].nt10 += 1
							if event.TLCat is 2 : res[ch_str]['ee'].nt01 += 1
							if event.TLCat is 3 : res[ch_str]['ee'].nt0  += 1

					if write_ResTree :
						# store oberservations in results tree
						self.ResTree_ObsPred[0] = 0
						self.ResTree_Weight [0] = 1.
						self.results_tree.Fill()

						# store fake predictions in results tree
						self.ResTree_ObsPred[0] = 1
						self.ResTree_Weight [0] = npf+nfp+nff
						self.results_tree.Fill()

				if write_ResTree and event.Flavor > 2 :
					self.ResTree_ObsPred[0] = 2

					# EM OS
					if event.Flavor is 4 :
						if event.TLCat is 0 : self.ResTree_Weight[0] = fbb
						if event.TLCat is 1 : self.ResTree_Weight[0] = fee

					# EE OS
					if event.Flavor is 5 :
						if event.TLCat is 0                     : self.ResTree_Weight[0] = 2*fbb
						if event.TLCat is 1 or event.TLCat is 2 : self.ResTree_Weight[0] = 2*feb
						if event.TLCat is 3                     : self.ResTree_Weight[0] = 2*fee

					self.results_tree.Fill()

			# GET RARE MC EVENTS
			if event.SType is 15 and event.TLCat == 0 and event.Flavor < 3 :
				if event.SName == 'WWTo2L2Nu' : continue # TODO: why?
				scale = event.PUWeight * event.HLTSF * self.lumi / self.samples[str(event.SName)].getLumi()

				if str(event.SName) not in rares :
					rares      [str(event.SName)] = {}
					rares_npass[str(event.SName)] = {}
					for ch_str, charge in self.charges.iteritems() :
						rares      [str(event.SName)][ch_str] = {}
						rares_npass[str(event.SName)][ch_str] = {}
						for chan in res[ch_str] :
							rares      [str(event.SName)][ch_str][chan] = 0.
							rares_npass[str(event.SName)][ch_str][chan] = 0

				for ch_str, charge in self.charges.iteritems() :
					if charge != 0 and event.Charge != charge : continue

					rares      [str(event.SName)][ch_str]['al'] += scale
					rares_npass[str(event.SName)][ch_str]['al'] += 1

					if event.Flavor is 0 :
						rares      [str(event.SName)][ch_str]['mm'] += scale
						rares_npass[str(event.SName)][ch_str]['mm'] += 1

					if event.Flavor is 1 :
						rares      [str(event.SName)][ch_str]['em'] += scale
						rares_npass[str(event.SName)][ch_str]['em'] += 1

					if event.Flavor is 2 :
						rares      [str(event.SName)][ch_str]['ee'] += scale
						rares_npass[str(event.SName)][ch_str]['ee'] += 1

				if write_ResTree :
					if   event.SName == 'WZTo3LNu' : self.ResTree_ObsPred[0] = 3
					elif event.SName == 'TTbarW'   : self.ResTree_ObsPred[0] = 4
					elif event.SName == 'TTbarZ'   : self.ResTree_ObsPred[0] = 5
					else                           : self.ResTree_ObsPred[0] = 6
					self.ResTree_Weight[0] = scale
					self.results_tree.Fill()

		# END LOOP OVER TREE

		if write_ResTree :
			self.results_file.Write()
			self.results_file.Close()

		sigfile.Close()

		# print number of events per sample in skimmed SigEvents tree
#		for name, n in nevents.iteritems() :
#			print '%12s: %8d' % (name, n)
#		for name, value in rares_mm.iteritems() :

		################
		# Observations #
		################

		print '[status] storing observations..'

		for ch_str in res :
			for chan in res[ch_str] :
				res[ch_str][chan].set_observations()

		#####################
		# Fakes predictions #
		#####################

		print '[status] calculating fake predictions..'

		FR.setNToyMCs(100)
		FR.setAddESyst(self.FakeESyst)

		# numbers from SSDLPlotter.cc just to check
#		mufratio_data = 0.040942; mufratio_data_e = 0.002156;
#		mupratio_data = 0.804292; mupratio_data_e = 0.001193;
#		elfratio_data = 0.069959; elfratio_data_e = 0.001419;
#		elpratio_data = 0.750609; elpratio_data_e = 0.001473;
#		FR.setMFRatio(mufratio_data, mufratio_data_e) # set error to pure statistical of ratio
#		FR.setEFRatio(elfratio_data, elfratio_data_e)
#		FR.setMPRatio(mupratio_data, mupratio_data_e)
#		FR.setEPRatio(elpratio_data, elpratio_data_e)

		# TODO: use these (ratios with ewk subtraction)
		FR.setMFRatio(self.fpr.MufRatio, self.fpr.MufRatioE) # set error to pure statistical of ratio
		FR.setEFRatio(self.fpr.ElfRatio, self.fpr.ElfRatioE)
		FR.setMPRatio(self.fpr.MupRatio, self.fpr.MupRatioE)
		FR.setEPRatio(self.fpr.ElpRatio, self.fpr.ElpRatioE)

		for ch_str in self.charges :

			FR.setMMNtl(res[ch_str]['mm'].nt2, res[ch_str]['mm'].nt10, res[ch_str]['mm'].nt01, res[ch_str]['mm'].nt0)
			FR.setEMNtl(res[ch_str]['em'].nt2, res[ch_str]['em'].nt10, res[ch_str]['em'].nt01, res[ch_str]['em'].nt0)
			FR.setEENtl(res[ch_str]['ee'].nt2, res[ch_str]['ee'].nt10, res[ch_str]['ee'].nt01, res[ch_str]['ee'].nt0)

			# store stat and syst errors
			res[ch_str]['mm'].npp_staterr = FR.getMMNppEStat(); res[ch_str]['mm'].npp_systerr = self.FakeESyst*res[ch_str]['mm'].npp;
			res[ch_str]['em'].npp_staterr = FR.getEMNppEStat(); res[ch_str]['em'].npp_systerr = self.FakeESyst*res[ch_str]['em'].npp;
			res[ch_str]['ee'].npp_staterr = FR.getEENppEStat(); res[ch_str]['ee'].npp_systerr = self.FakeESyst*res[ch_str]['ee'].npp;

			res[ch_str]['mm'].npf_staterr = FR.getMMNpfEStat(); res[ch_str]['mm'].npf_systerr = self.FakeESyst*res[ch_str]['mm'].npf;
			res[ch_str]['em'].npf_staterr = FR.getEMNpfEStat(); res[ch_str]['em'].npf_systerr = self.FakeESyst*res[ch_str]['em'].npf;
			res[ch_str]['ee'].npf_staterr = FR.getEENpfEStat(); res[ch_str]['ee'].npf_systerr = self.FakeESyst*res[ch_str]['ee'].npf;

			res[ch_str]['mm'].nfp_staterr = FR.getMMNfpEStat(); res[ch_str]['mm'].nfp_systerr = self.FakeESyst*res[ch_str]['mm'].nfp;
			res[ch_str]['em'].nfp_staterr = FR.getEMNfpEStat(); res[ch_str]['em'].nfp_systerr = self.FakeESyst*res[ch_str]['em'].nfp;
			res[ch_str]['ee'].nfp_staterr = FR.getEENfpEStat(); res[ch_str]['ee'].nfp_systerr = self.FakeESyst*res[ch_str]['ee'].nfp;

			res[ch_str]['mm'].nff_staterr = FR.getMMNffEStat(); res[ch_str]['mm'].nff_systerr = self.FakeESyst*res[ch_str]['mm'].nff;
			res[ch_str]['em'].nff_staterr = FR.getEMNffEStat(); res[ch_str]['em'].nff_systerr = self.FakeESyst*res[ch_str]['em'].nff;
			res[ch_str]['ee'].nff_staterr = FR.getEENffEStat(); res[ch_str]['ee'].nff_systerr = self.FakeESyst*res[ch_str]['ee'].nff;

			# store fake predictions
			res[ch_str]['mm'].set_fakePredictions()
			res[ch_str]['em'].set_fakePredictions()
			res[ch_str]['ee'].set_fakePredictions()
			res[ch_str]['al'].fake = res[ch_str]['mm'].fake + res[ch_str]['em'].fake + res[ch_str]['ee'].fake
			res[ch_str]['al'].fake_err = math.sqrt(FR.getTotEStat()  *FR.getTotEStat()   + self.FakeESyst2*res[ch_str]['al'].fake*res[ch_str]['al'].fake)
			res[ch_str]['mm'].fake_err = math.sqrt(FR.getMMTotEStat()*FR.getMMTotEStat() + self.FakeESyst2*res[ch_str]['mm'].fake*res[ch_str]['mm'].fake)
			res[ch_str]['em'].fake_err = math.sqrt(FR.getEMTotEStat()*FR.getEMTotEStat() + self.FakeESyst2*res[ch_str]['em'].fake*res[ch_str]['em'].fake)
			res[ch_str]['ee'].fake_err = math.sqrt(FR.getEETotEStat()*FR.getEETotEStat() + self.FakeESyst2*res[ch_str]['ee'].fake*res[ch_str]['ee'].fake)

		#########
		# Rares #
		#########

		print '[status] storing all rare SM processes..'

		for ch_str in self.charges :

			for chan in res[ch_str] :

				wz   = 0.; wz_staterr2   = 0.;
				ttw  = 0.; ttw_staterr2  = 0.;
				ttz  = 0.; ttz_staterr2  = 0.;
				rare = 0.; rare_staterr2 = 0.;

				for s in rares :
					scale = self.lumi / self.samples[s].getLumi()
					#scale = HLTSF * self.lumi / self.samples[s].getLumi() # this reproduces the bug in the SSDLPlotter.cc
					staterr = scale * self.samples[s].getError(rares_npass[s][ch_str][chan])

					if s == 'WZTo3LNu' :
						wz += rares[s][ch_str][chan]
						wz_staterr2 += staterr * staterr

					elif s == 'TTbarW' :
						ttw += rares[s][ch_str][chan]
						ttw_staterr2 += staterr * staterr

					elif s == 'TTbarZ' :
						ttz += rares[s][ch_str][chan]
						ttz_staterr2 += staterr * staterr

					else :
						rare += rares[s][ch_str][chan]
						rare_staterr2 += staterr * staterr

						res[ch_str][chan].rares[s] = rares[s][ch_str][chan]
						res[ch_str][chan].rares_staterr[s] = staterr

				# store WZ yields
				res[ch_str][chan].wz      = wz
				res[ch_str][chan].wz_err  = math.sqrt(wz_staterr2 + self.WZESyst2 * wz * wz)

				# store ttW mc yields
				res[ch_str][chan].ttw     = ttw
				res[ch_str][chan].ttw_err = math.sqrt(ttw_staterr2 + self.TTWESyst2 * ttw * ttw)

				# store ttZ mc yields
				res[ch_str][chan].ttz     = ttz
				res[ch_str][chan].ttz_err = math.sqrt(ttz_staterr2 + self.TTZESyst2 * ttz * ttz)

				# store rare mc yields
				res[ch_str][chan].rare     = rare
				res[ch_str][chan].rare_err = math.sqrt(rare_staterr2 + self.RareESyst2 * rare * rare)

				# store ttW/Z mc yields
				res[ch_str][chan].set_ttwzPredictions()

		#####################
		# Total backgrounds #
		#####################

		for ch_str in res :
			for chan in res[ch_str] :
				res[ch_str][chan].set_totBackground()

		##########################################
		# print all observations and predictions #
		##########################################

#		yields.printout()
		self.print_results(res['al'])
		self.print_results(res['++'])

		return res


	def print_results(self, res) :
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
				res['mm'].rares[s], res['mm'].rares_staterr[s], self.RareESyst*res['mm'].rares[s],
				res['em'].rares[s], res['em'].rares_staterr[s], self.RareESyst*res['em'].rares[s],
				res['ee'].rares[s], res['ee'].rares_staterr[s], self.RareESyst*res['ee'].rares[s],
				res['al'].rares[s], res['al'].rares_staterr[s], self.RareESyst*res['al'].rares[s])


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
		print "----------------------------------------------------------------------------------------------"


	def make_datacard(self, results, chan, charge) :
		'''
		takes a nested dictionary of result objects as input:
		results[SYSTFLAG][CHARGE][FLAVOR]

		writes datacard in the combine tool format
		'''

		datacard_name = 'datacard_ssdl_ttW_' + results['Normal'][charge][chan].chan_str + '.txt'
		datacard_path = self.path + 'datacards/'
		helper.mkdir(datacard_path)
		print '[status] writing %s' % datacard_name
		with open(datacard_path + datacard_name, 'w') as file :
			timestamp = time.asctime()
			file.write('#=========================================================================================\n')
			file.write('# Systematics table for ttW analysis, same-sign channel, subchannels\n')
			file.write('# Generated on: %s\n' % (timestamp))
			file.write('# Copy between the dashed lines for datacard\n')
			file.write('#-----------------------------------------------------------------------------------------\n')
			file.write('imax 1\n')
			file.write('jmax 5\n')
			file.write('kmax *\n')
			file.write('\n')
			file.write('bin\t\t%s\n' % (results['Normal'][charge][chan].chan_str))
			##	if (gFullDataBlind)
			##		fOUTSTREAM << Form("observation\t%d\t%d\t%d\t%d\t%d\t%d", 999, 999, 999, 999, 999, 999) << endl;
			##	else
			file.write('observation\t%d\n' % (results['Normal'][charge][chan].obs))
			file.write('\n')
			file.write('bin\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\n' % (
				results['Normal'][charge][chan].chan_str,
				results['Normal'][charge][chan].chan_str,
				results['Normal'][charge][chan].chan_str,
				results['Normal'][charge][chan].chan_str,
				results['Normal'][charge][chan].chan_str,
				results['Normal'][charge][chan].chan_str))
			file.write('process\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare\n')
			file.write('process\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5\n')
			file.write('rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\n' % (
				results['Normal'][charge][chan].ttw,
				results['Normal'][charge][chan].ttz,
				results['Normal'][charge][chan].fake,
				results['Normal'][charge][chan].cmid,
				results['Normal'][charge][chan].wz,
				results['Normal'][charge][chan].rare))
			file.write('\n')
			file.write('#syst\n')

			file.write('lumi     lnN\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\n' % (
				self.lumi_syst,
				self.lumi_syst,
				self.lumi_syst,
				self.lumi_syst))

			file.write(    'bgUncttz lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-\n' % (1. + results['Normal'][charge][chan].ttz_err  / results['Normal'][charge][chan].ttz ))
			file.write(    'bgUncfak lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-\n' % (1. + results['Normal'][charge][chan].fake_err / results['Normal'][charge][chan].fake))
			if chan != 'mm' :
				file.write('bgUnccmi lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-\n' % (1. + results['Normal'][charge][chan].cmid_err / results['Normal'][charge][chan].cmid))
			file.write(    'bgUncwz  lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\n' % (1. + results['Normal'][charge][chan].wz_err   / results['Normal'][charge][chan].wz  ))
			file.write(    'bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f\n' % (1. + results['Normal'][charge][chan].rare_err / results['Normal'][charge][chan].rare))

			if 'LepUp' in results and 'LepDown' in results :
				file.write('lept     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\n' % (
					results['LepDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['LepUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['LepDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['LepUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['LepDown'][charge][chan].wz   / results['Normal'][charge][chan].wz  ,
					results['LepUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz  ,
					results['LepDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
					results['LepUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare))
			else :
				print '[WARNING] LepUp/Down systematic not found!'

			if 'MuUp' in results and 'MuDown' in results :
				file.write('muon     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\n' % (
					results['MuDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['MuUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['MuDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['MuUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['MuDown'][charge][chan].wz   / results['Normal'][charge][chan].wz  ,
					results['MuUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz  ,
					results['MuDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
					results['MuUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare))
			else :
				print '[WARNING] MuUp/Down systematic not found!'

			if 'ElUp' in results and 'ElDown' in results :
				file.write('elec     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\n' % (
					results['ElDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['ElUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['ElDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['ElUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['ElDown'][charge][chan].wz   / results['Normal'][charge][chan].wz  ,
					results['ElUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz  ,
					results['ElDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
					results['ElUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare))
			else :
				print '[WARNING] ElUp/Down systematic not found!'

			file.write('leptrig  lnN\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\n' % (
				self.ltrig_syst,
				self.ltrig_syst,
				self.ltrig_syst,
				self.ltrig_syst))

			if 'BUp' in results and 'BDown' in results :
				file.write('btag     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\n' % (
					results['BDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['BUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['BDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['BUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['BDown'][charge][chan].wz   / results['Normal'][charge][chan].wz,
					results['BUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz,
					results['BDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
					results['BUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare))
			else :
				print '[WARNING] BUp/Down systematic not found!'

			if 'JetUp' in results and 'JetDown' in results :
				file.write('jes      lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\n' % (
					results['JetDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['JetUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['JetDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['JetUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['JetDown'][charge][chan].wz   / results['Normal'][charge][chan].wz,
					results['JetUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz,
					results['JetDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
					results['JetUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare))
			else :
				print '[WARNING] JetUp/Down systematic not found!'

			if 'JetSmearUp' in results and 'JetSmearDown' in results :
				file.write('jer      lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\n' % (
					results['JetSmearDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['JetSmearUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['JetSmearDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['JetSmearUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['JetSmearDown'][charge][chan].wz   / results['Normal'][charge][chan].wz,
					results['JetSmearUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz,
					results['JetSmearDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
					results['JetSmearUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare))
			else :
				print '[WARNING] JetSmearUp/Down systematic not found!'

			if 'PileupUp' in results and 'PileupDown' in results :
				file.write('pu       lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\n' % (
					results['PileupDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['PileupUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['PileupDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['PileupUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
					results['PileupDown'][charge][chan].wz   / results['Normal'][charge][chan].wz,
					results['PileupUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz,
					results['PileupDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
					results['PileupUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare))
			else :
				print '[WARNING] PileupUp/Down systematic not found!'

			##//	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
			##//					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn,
			##//					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn,
			##//					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn,
			##//					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
			file.write('scale    lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\n' % (
				self.scale_syst_up, self.scale_syst_dn, self.scale_syst_up, self.scale_syst_dn, self.scale_syst_up, self.scale_syst_dn, self.scale_syst_up, self.scale_syst_dn))

			file.write('tmass    lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t-\t\t-\n' % (
				self.tmass_syst_up, self.tmass_syst_dn, self.tmass_syst_up, self.tmass_syst_dn))

			file.write('gen      lnN\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\n' % (self.gen_syst))

			file.write('pdf      lnN\t%5.3f\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\n' % (self.pdf_syst, self.wz_pdf_syst))


	def plot_predictions(self, treepath, sel) :
		vars = []
		vars.append('HT'    )
		vars.append('MET'   )
		vars.append('NJ'    )
		vars.append('NbJmed')
		vars.append('pT1'   )
		vars.append('pT2'   )
		vars.append('Mll'   )
		vars.append('NVrtx' )
		vars.append('minMT' )
		vars.append('M3'    )
		vars.append('Int'   )
		vars.append('CFChan')  # ChMisID prediction is not totally correct and it's statistical uncertainty a bit too large. Since not all OS events are considered and diveded by 2, but only the ones in the according charge channel.
		vars.append('Charge')  # ChMisID prediction is not totally correct and it's statistical uncertainty a bit too large. Since not all OS events are considered and diveded by 2, but only the ones in the according charge channel.

		print '[status] open Results tree from %s' % (treepath)
		resfile = ROOT.TFile.Open(treepath, 'READ')
		restree = resfile.Get('Results')

		print '[status] getting plots from %s tree for' % (restree.GetName()), vars

		for var in vars :
			histo_bins = config.get_histoBins(var, sel)
			self.plot_ObsPred(restree, sel, var, histo_bins)

		resfile.Close()


	def plot_ObsPred(self, tree, sel, var, settings) :

		nbins = settings['nbins']
		min   = settings['min'  ]
		max   = settings['max'  ]

		histos = {}

		chargeFactor = 1.
		if sel.charge != 0 : chargeFactor = 0.5

		# loop over results tree to get names of rare samples
		rares = []
		wz_str  = ''
		ttw_str = ''
		ttz_str = ''
		for event in tree :
			if event.ObsPred == 6 and str(event.SName) not in rares :
				rares.append(str(event.SName))
			if event.ObsPred == 3 and wz_str  == '' : wz_str  = str(event.SName)
			if event.ObsPred == 4 and ttw_str == '' : ttw_str = str(event.SName)
			if event.ObsPred == 5 and ttz_str == '' : ttz_str = str(event.SName)

		############################
		# SETUP AND GET HISTOGRAMS #
		############################

		if   var == 'pT1'   : var_str = 'TMath::Max(pT1,pT2)'
		elif var == 'pT2'   : var_str = 'TMath::Min(pT1,pT2)'
		elif var == 'minMT' : var_str = 'TMath::Min(MTLep1,MTLep2)'
		elif var == 'Int'   : var_str = 'Flavor'
		else                : var_str = var

		# get histograms with nPass for rares, wz, ttz and ttw
		h_rares_npass = {}
		for rare in rares :
			h_rare_name = 'h_'+rare+'_npass_'+var; h_rares_npass[rare] = ROOT.TH1D(h_rare_name, h_rare_name, nbins, min, max)
			tree.Draw(var_str+'>>'+h_rare_name, '(ObsPred == 6 && SName == \"%s\" && %s)' % (rare, sel.get_selectionString(ResTree = True)), 'goff')
		h_wz_name    = 'h_wz_npass'    + var + sel.name; h_wz_npass  = ROOT.TH1D(h_wz_name   , h_wz_name   , nbins, min, max)
		h_ttz_name   = 'h_ttz_npass'   + var + sel.name; h_ttz_npass = ROOT.TH1D(h_ttz_name  , h_ttz_name  , nbins, min, max)
		h_ttw_name   = 'h_ttw_npass'   + var + sel.name; h_ttw_npass = ROOT.TH1D(h_ttw_name  , h_ttw_name  , nbins, min, max)
		tree.Draw(var_str+'>>'+h_wz_name   , '(ObsPred == 3 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_ttz_name  , '(ObsPred == 5 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_ttw_name  , '(ObsPred == 4 && %s)' % sel.get_selectionString(ResTree = True), 'goff')

		# setup histograms for observation, fakes, chmid, wz, ttw, ttz and rare
		h_obs_name   = 'h_obs_'   + var + sel.name; histos['obs'  ] = ROOT.TH1D(h_obs_name  , h_obs_name  , nbins, min, max)
		h_fake_name  = 'h_fake_'  + var + sel.name; histos['fake' ] = ROOT.TH1D(h_fake_name , h_fake_name , nbins, min, max)
		h_chmid_name = 'h_chmid_' + var + sel.name; histos['chmid'] = ROOT.TH1D(h_chmid_name, h_chmid_name, nbins, min, max)
		h_rare_name  = 'h_rare_'  + var + sel.name; histos['rare' ] = ROOT.TH1D(h_rare_name , h_rare_name , nbins, min, max)
		h_wz_name    = 'h_wz_'    + var + sel.name; histos['wz'   ] = ROOT.TH1D(h_wz_name   , h_wz_name   , nbins, min, max)
		h_ttz_name   = 'h_ttz_'   + var + sel.name; histos['ttz'  ] = ROOT.TH1D(h_ttz_name  , h_ttz_name  , nbins, min, max)
		h_ttw_name   = 'h_ttw_'   + var + sel.name; histos['ttw'  ] = ROOT.TH1D(h_ttw_name  , h_ttw_name  , nbins, min, max)
		h_bgtot_name = 'h_bgtot_' + var + sel.name; histos['bgtot'] = ROOT.TH1D(h_bgtot_name, h_bgtot_name, nbins, min, max)
		h_pred_name  = 'h_pred_'  + var + sel.name; histos['pred' ] = ROOT.TH1D(h_pred_name , h_pred_name , nbins, min, max)
		h_mm_nt2_name  = 'h_mm_nt2_'  + var + sel.name; h_mm_nt2_npass  = ROOT.TH1D(h_mm_nt2_name , h_mm_nt2_name , nbins, min, max)
		h_mm_nt10_name = 'h_mm_nt10_' + var + sel.name; h_mm_nt10_npass = ROOT.TH1D(h_mm_nt10_name, h_mm_nt10_name, nbins, min, max)
		h_mm_nt01_name = 'h_mm_nt01_' + var + sel.name; h_mm_nt01_npass = ROOT.TH1D(h_mm_nt01_name, h_mm_nt01_name, nbins, min, max)
		h_mm_nt0_name  = 'h_mm_nt0_'  + var + sel.name; h_mm_nt0_npass  = ROOT.TH1D(h_mm_nt0_name , h_mm_nt0_name , nbins, min, max)
		h_em_nt2_name  = 'h_em_nt2_'  + var + sel.name; h_em_nt2_npass  = ROOT.TH1D(h_em_nt2_name , h_em_nt2_name , nbins, min, max)
		h_em_nt10_name = 'h_em_nt10_' + var + sel.name; h_em_nt10_npass = ROOT.TH1D(h_em_nt10_name, h_em_nt10_name, nbins, min, max)
		h_em_nt01_name = 'h_em_nt01_' + var + sel.name; h_em_nt01_npass = ROOT.TH1D(h_em_nt01_name, h_em_nt01_name, nbins, min, max)
		h_em_nt0_name  = 'h_em_nt0_'  + var + sel.name; h_em_nt0_npass  = ROOT.TH1D(h_em_nt0_name , h_em_nt0_name , nbins, min, max)
		h_ee_nt2_name  = 'h_ee_nt2_'  + var + sel.name; h_ee_nt2_npass  = ROOT.TH1D(h_ee_nt2_name , h_ee_nt2_name , nbins, min, max)
		h_ee_nt10_name = 'h_ee_nt10_' + var + sel.name; h_ee_nt10_npass = ROOT.TH1D(h_ee_nt10_name, h_ee_nt10_name, nbins, min, max)
		h_ee_nt01_name = 'h_ee_nt01_' + var + sel.name; h_ee_nt01_npass = ROOT.TH1D(h_ee_nt01_name, h_ee_nt01_name, nbins, min, max)
		h_ee_nt0_name  = 'h_ee_nt0_'  + var + sel.name; h_ee_nt0_npass  = ROOT.TH1D(h_ee_nt0_name , h_ee_nt0_name , nbins, min, max)
		h_nt2_em_BB_os_name = 'h_nt2_em_BB_os_' + var + sel.name; h_nt2_em_BB_os = ROOT.TH1D(h_nt2_em_BB_os_name, h_nt2_em_BB_os_name, nbins, min, max)
		h_nt2_em_EE_os_name = 'h_nt2_em_EE_os_' + var + sel.name; h_nt2_em_EE_os = ROOT.TH1D(h_nt2_em_EE_os_name, h_nt2_em_EE_os_name, nbins, min, max)
		h_nt2_ee_BB_os_name = 'h_nt2_ee_BB_os_' + var + sel.name; h_nt2_ee_BB_os = ROOT.TH1D(h_nt2_ee_BB_os_name, h_nt2_ee_BB_os_name, nbins, min, max)
		h_nt2_ee_EB_os_name = 'h_nt2_ee_EB_os_' + var + sel.name; h_nt2_ee_EB_os = ROOT.TH1D(h_nt2_ee_EB_os_name, h_nt2_ee_EB_os_name, nbins, min, max)
		h_nt2_ee_EE_os_name = 'h_nt2_ee_EE_os_' + var + sel.name; h_nt2_ee_EE_os = ROOT.TH1D(h_nt2_ee_EE_os_name, h_nt2_ee_EE_os_name, nbins, min, max)

		# getting histograms from results tree
		tree.Draw(var_str+'>>'+h_obs_name  ,    'Weight*(ObsPred == 0 && TLCat == 0 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_fake_name ,    'Weight*(ObsPred == 1 && %s)'               % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_chmid_name, '%f*Weight*(ObsPred == 2 && %s)' %(chargeFactor, sel.get_selectionString(ResTree = True, OS_data = (0,-1))), 'goff')
		tree.Draw(var_str+'>>'+h_rare_name ,    'Weight*(ObsPred == 6 && TLCat == 0 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_wz_name   ,    'Weight*(ObsPred == 3 && TLCat == 0 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_ttz_name  ,    'Weight*(ObsPred == 5 && TLCat == 0 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_ttw_name  ,    'Weight*(ObsPred == 4 && TLCat == 0 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_mm_nt2_name ,         '(ObsPred == 0 && TLCat == 0 && Flavor == 0 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_mm_nt10_name,         '(ObsPred == 0 && TLCat == 1 && Flavor == 0 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_mm_nt01_name,         '(ObsPred == 0 && TLCat == 2 && Flavor == 0 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_mm_nt0_name ,         '(ObsPred == 0 && TLCat == 3 && Flavor == 0 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_em_nt2_name ,         '(ObsPred == 0 && TLCat == 0 && Flavor == 1 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_em_nt10_name,         '(ObsPred == 0 && TLCat == 1 && Flavor == 1 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_em_nt01_name,         '(ObsPred == 0 && TLCat == 2 && Flavor == 1 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_em_nt0_name ,         '(ObsPred == 0 && TLCat == 3 && Flavor == 1 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_ee_nt2_name ,         '(ObsPred == 0 && TLCat == 0 && Flavor == 2 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_ee_nt10_name,         '(ObsPred == 0 && TLCat == 1 && Flavor == 2 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_ee_nt01_name,         '(ObsPred == 0 && TLCat == 2 && Flavor == 2 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_ee_nt0_name ,         '(ObsPred == 0 && TLCat == 3 && Flavor == 2 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_nt2_em_BB_os_name,    '(ObsPred == 2 && TLCat == 0 && Flavor == 1 && %s)' % sel.get_selectionString(ResTree = True, OS_data = (0,-1)), 'goff')
		tree.Draw(var_str+'>>'+h_nt2_em_EE_os_name,    '(ObsPred == 2 && TLCat == 1 && Flavor == 1 && %s)' % sel.get_selectionString(ResTree = True, OS_data = (0,-1)), 'goff')
		tree.Draw(var_str+'>>'+h_nt2_ee_BB_os_name,    '(ObsPred == 2 && TLCat == 0 && Flavor == 2 && %s)' % sel.get_selectionString(ResTree = True, OS_data = (0,-1)), 'goff')
		tree.Draw(var_str+'>>'+h_nt2_ee_EB_os_name,    '(ObsPred == 2 && (TLCat == 1 || TLCat == 2) && Flavor == 2 && %s)' % sel.get_selectionString(ResTree = True, OS_data = (0,-1)), 'goff')
		tree.Draw(var_str+'>>'+h_nt2_ee_EE_os_name,    '(ObsPred == 2 && TLCat == 3 && Flavor == 2 && %s)' % sel.get_selectionString(ResTree = True, OS_data = (0,-1)), 'goff')

		# adding background predictions
		histos['bgtot'].Add(histos['fake' ])
		histos['bgtot'].Add(histos['chmid'])
		histos['bgtot'].Add(histos['rare' ])
		histos['bgtot'].Add(histos['wz'   ])
		histos['bgtot'].Add(histos['ttz'  ])

		# adding predictions
		histos['pred'].Add(histos['fake' ])
		histos['pred'].Add(histos['chmid'])
		histos['pred'].Add(histos['rare' ])
		histos['pred'].Add(histos['wz'   ])
		histos['pred'].Add(histos['ttz'  ])
		histos['pred'].Add(histos['ttw'  ])

		##########
		# ERRORS #
		##########

		FR = ROOT.FakeRatios()
		FR.setNToyMCs(100)
		FR.setAddESyst(self.FakeESyst)

		FR.setMFRatio(self.fpr.MufRatio, self.fpr.MufRatioE) # set error to pure statistical of ratio
		FR.setEFRatio(self.fpr.ElfRatio, self.fpr.ElfRatioE)
		FR.setMPRatio(self.fpr.MupRatio, self.fpr.MupRatioE)
		FR.setEPRatio(self.fpr.ElpRatio, self.fpr.ElpRatioE)

		(fbb, fbbE) = self.calculateChMisIdProb(self.get_samples('DoubleEle'), 'BB', self.chmid_sf)
		(feb, febE) = self.calculateChMisIdProb(self.get_samples('DoubleEle'), 'EB', self.chmid_sf)
		(fee, feeE) = self.calculateChMisIdProb(self.get_samples('DoubleEle'), 'EE', self.chmid_sf)

		for bin in range(1, histos['bgtot'].GetNbinsX()+1) :
			err2 = 0.

			# fakes
			FR.setMMNtl(h_mm_nt2_npass.GetBinContent(bin), h_mm_nt10_npass.GetBinContent(bin), h_mm_nt01_npass.GetBinContent(bin), h_mm_nt0_npass.GetBinContent(bin))
			FR.setEMNtl(h_em_nt2_npass.GetBinContent(bin), h_em_nt10_npass.GetBinContent(bin), h_em_nt01_npass.GetBinContent(bin), h_em_nt0_npass.GetBinContent(bin))
			FR.setEENtl(h_ee_nt2_npass.GetBinContent(bin), h_ee_nt10_npass.GetBinContent(bin), h_ee_nt01_npass.GetBinContent(bin), h_ee_nt0_npass.GetBinContent(bin))

			fake_nPass = histos['fake'].GetBinContent(bin)
			fake_syst2 = self.FakeESyst2 * fake_nPass*fake_nPass
			fake_stat2 = FR.getTotEStat() * FR.getTotEStat()
			err2 += fake_syst2 + fake_stat2
#			print 'bin %d: %f +/- %f (syst) +/- %f (stat)' % (bin, fake_nPass, math.sqrt(fake_syst2), math.sqrt(fake_stat2))

			# charge mis ID
			nt2_em_BB_os = h_nt2_em_BB_os.GetBinContent(bin)
			nt2_em_EE_os = h_nt2_em_EE_os.GetBinContent(bin)
			nt2_ee_BB_os = h_nt2_ee_BB_os.GetBinContent(bin)
			nt2_ee_EB_os = h_nt2_ee_EB_os.GetBinContent(bin)
			nt2_ee_EE_os = h_nt2_ee_EE_os.GetBinContent(bin)
				# Simple error propagation assuming error on number of events is sqrt(N)
			nt2_ee_chmid    = chargeFactor * (2*fbb* nt2_ee_BB_os                           + 2*fee*nt2_ee_EE_os                      + 2*feb*nt2_ee_EB_os)
			nt2_ee_chmid_e1 = chargeFactor * math.sqrt( 4*fbb*fbb * FR.getEStat2(nt2_ee_BB_os)  + 4*fee*fee * FR.getEStat2(nt2_ee_EE_os)    + 4*feb*feb * FR.getEStat2(nt2_ee_EB_os) )  # stat only
			nt2_ee_chmid_e2 = chargeFactor * math.sqrt( 4*fbbE*fbbE * nt2_ee_BB_os*nt2_ee_BB_os + 4*feeE*feeE * nt2_ee_EE_os*nt2_ee_EE_os   + 4*febE*febE * nt2_ee_EB_os*nt2_ee_EB_os + self.ChMisESyst2 * nt2_ee_chmid*nt2_ee_chmid/(chargeFactor*chargeFactor) )  # syst only

			nt2_em_chmid    = chargeFactor * (fbb * nt2_em_BB_os + fee * nt2_em_EE_os)
			nt2_em_chmid_e1 = chargeFactor * math.sqrt( fbb*fbb * FR.getEStat2(nt2_em_BB_os)  + fee*fee*FR.getEStat2(nt2_em_EE_os) )
			nt2_em_chmid_e2 = chargeFactor * math.sqrt( fbbE*fbbE * nt2_em_BB_os*nt2_em_BB_os + feeE*feeE * nt2_em_EE_os*nt2_em_EE_os + 0.25*self.ChMisESyst2*nt2_em_chmid*nt2_em_chmid/(chargeFactor*chargeFactor) )
			err2 += nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2
#			print bin, histos['chmid'].GetBinContent(bin), nt2_em_chmid + nt2_ee_chmid

			# wz
			scale = self.lumi / self.samples[wz_str].getLumi()
			wz_nPass = h_wz_npass.GetBinContent(bin)
			wz_stat2 = scale*scale * self.samples[wz_str].getError2(wz_nPass)
			wz_yield = histos['wz'].GetBinContent(bin)
			wz_syst2 = self.WZESyst2 * wz_yield*wz_yield
			err2 += wz_syst2 + wz_stat2

			# ttz
			scale = self.lumi / self.samples[ttz_str].getLumi()
			ttz_nPass = h_ttz_npass.GetBinContent(bin)
			ttz_stat2 = scale*scale * self.samples[ttz_str].getError2(ttz_nPass)
			ttz_yield = histos['ttz'].GetBinContent(bin)
			ttz_syst2 = self.TTZESyst2 * ttz_yield*ttz_yield
			err2 += ttz_syst2 + ttz_stat2

#			# ttz
#			scale = self.lumi / self.samples[ttw_str].getLumi()
#			ttw_nPass = h_ttw_npass.GetBinContent(bin)
#			ttw_stat2 = scale*scale * self.samples[ttw_str].getError2(ttw_nPass)
#			ttw_yield = histos['ttw'].GetBinContent(bin)
#			ttw_syst2 = self.ttwESyst2 * ttw_yield*ttw_yield
#			err2 += ttw_syst2 + ttw_stat2

			# rare
			rare_yield = histos['rare'].GetBinContent(bin)
			rare_syst2 = self.RareESyst2 * rare_yield*rare_yield
			rare_stat2 = 0.
			for rare, h_rare_npass in h_rares_npass.iteritems() :
				scale = self.lumi / self.samples[rare].getLumi()
				rare_nPass = h_rare_npass.GetBinContent(bin)
				rare_stat2 += scale*scale * self.samples[rare].getError2(rare_nPass)
			err2 += rare_syst2 + rare_stat2

			histos['pred'].SetBinError(bin, math.sqrt(err2))

		histofile = ROOT.TFile.Open(self.path + 'SSDLHistos.root', 'RECREATE')
		histofile.cd()
		for process, histo in histos.iteritems() :
			histo.Write()
		histofile.Close()

		cms_label = 2
		if sel.flavor == -2 : prefix = 'OS'
		if sel.flavor == -1 : prefix = 'ALL'
		if sel.flavor ==  0 : prefix = 'MM'
		if sel.flavor ==  1 : prefix = 'EM'
		if sel.flavor ==  2 : prefix = 'EE'
		charge_str = ''
		if sel.charge != 0 :
			if sel.charge > 0 : charge_str = '^{+}'
			if sel.charge < 0 : charge_str = '^{-}'
		if sel.name == 'final' : cms_label = 0
		cms_label = 0
		pl = ttvplot.ttvplot(self.path + 'ObsPredPlots/%s/'%sel.name, '2L', self.lumi, cms_label)
		pl.save_plot(histos, var, sel.name, prefix, charge_str)
#		raw_input('ok? ')
#		raw_input('ok? ')


	def make_KinPlots(self) :
		foo = 0


if __name__ == '__main__' :
	args = sys.argv

	if ('--help' in args) or ('-h' in args) or ('-d' not in args) or ('-c' not in args) :
		print 'usage: ..'
		sys.exit(1)

	if ('-d' in args) and (args[args.index('-d')+1] != '') :
		path = str(args[args.index('-d')+1])
		print path

	if ('-c' in args) and (args[args.index('-c')+1] != '') :
		cardfile = str(args[args.index('-c')+1])
		print cardfile

	plotter = plotter(path)
	plotter.do_analysis(cardfile)
