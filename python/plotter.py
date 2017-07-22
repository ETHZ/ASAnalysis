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
import copy
import ttvStyle
import ConfigParser
import fakeratios


class plotter :
	'''the plotter reads sigtree and produces plots'''

	def __init__(self, path, cardfile, sel_path = '') :
		print '[status] initialize plotter..'
		self.path = path
		if not self.path.endswith('/') : self.path += '/'
		self.ssdlfile = ROOT.TFile.Open(self.path + 'SSDLYields.root', 'READ')
		if self.ssdlfile == None :
			sys.exit(1)

		# samples
		self.samples = self.readDatacard(cardfile)

		# selections
		if sel_path == '' :
			sel_path = '%s/selections.cfg' % os.path.dirname(os.path.realpath(__file__))
		self.selections = self.read_selections(sel_path)

		# ratios
		self.fpr = ratios.ratios(self.path, self.samples)

		# charge mis-ID
		self.chmid_sf = 1.

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
		self.systflags = {}
		self.systflags['Normal']       = 0
		self.systflags['JetUp']        = 1
		self.systflags['JetDown']      = 2
		self.systflags['JetSmear']     = 3
		self.systflags['BUp']          = 4
		self.systflags['BDown']        = 5
		self.systflags['LepUp']        = 6
		self.systflags['LepDown']      = 7
#		self.systflags['METUp']        = 8
#		self.systflags['METDown']      = 9
		self.systflags['PileupUp']     = 10
		self.systflags['PileupDown']   = 11
		self.systflags['JetSmearUp']   = 12
		self.systflags['JetSmearDown'] = 13
		self.systflags['MuUp']         = 14
		self.systflags['MuDown']       = 15
		self.systflags['ElUp']         = 16
		self.systflags['ElDown']       = 17

		# systematic uncertainties for datacard
		self.systematics = {}
		self.systematics['lumi'    ] = 1.026
		self.systematics['scale_up'] = 1.027
		self.systematics['scale_dn'] = 0.980
		self.systematics['tmass_up'] = 1.019
		self.systematics['tmass_dn'] = 0.983
		self.systematics['ltrig'   ] = 1.030
		self.systematics['pdf'     ] = 1.015
		self.systematics['gen'     ] = 1.050
		self.systematics['wz_pdf'  ] = 1.041

		# charge strings
		self.charges = {}
		for charge in range(1, -2, -1) :
			self.charges[self.get_chargeString(charge)] = charge

		# random variable
		self.rand = ROOT.TRandom3(0)


	def do_analysis(self, IntPred = True, DiffPred = False, IntMC = False, DiffMC = False, RatioPlots = False, RatioControlPlots = False, FakeClosure = False, SLPlots = False) :
		print '[status] starting analysis..'

		# make table of samples
		tables.make_SampleTable(self.path, [self.samples[name] for name in self.get_samples('Rare')], 'Irreducible')
		tables.make_SampleTable(self.path, [self.samples[name] for name in self.get_samples('MC'  )], 'MC'         )

		# charge mis-ID scale factor
		self.chmid_sf = 1.62

		# get fake and prompt ratios
#		EWK_SF = {}
#		EWK_SF['el']   = self.get_EWK_SF('el')
#		EWK_SF['mu17'] = self.get_EWK_SF('mu17')
#		EWK_SF['mu24'] = self.get_EWK_SF('mu24')
		if RatioPlots :
			self.fpr.fill_ratios(self.get_samples('SingleDoubleMu'), self.get_samples('DoubleEle'), 0, False)
			self.fpr.fill_ratios(self.get_samples('MC')            , self.get_samples('MC')       , 1, False)
			#self.fpr.fill_ratios(self.get_samples('QCD')            , self.get_samples('QCD')       , 1, False)
			self.fpr.plot_ratios('Ratios_uncorrected')

		self.fpr.fill_ratios(self.get_samples('SingleDoubleMu'), self.get_samples('DoubleEle'), 0, True)
		self.fpr.fill_ratios(self.get_samples('MC')            , self.get_samples('MC')       , 1, True)

		if RatioControlPlots != False :
			if RatioControlPlots == 'el' :
				self.fpr.make_controlPlots(self.get_samples('DoubleEle'), RatioControlPlots)
				self.fpr.make_controlPlots(self.get_samples('MC'       ), RatioControlPlots, datamc = 1)
			if RatioControlPlots == 'mu17' :
				self.fpr.make_controlPlots(self.get_samples('DoubleMu' ), RatioControlPlots)
			if RatioControlPlots == 'mu24' :
				self.fpr.make_controlPlots(self.get_samples('SingleMu' ), RatioControlPlots)

		if RatioPlots :
			self.fpr.plot_ratios()

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


		if DiffMC :
#			self.plot_DiffMC(self.selections['1J0bJ'    ])
#			self.plot_DiffMC(self.selections['2J0bJ'    ])
#			self.plot_DiffMC(self.selections['2J0bJOS'        ])
			self.plot_DiffMC(self.selections['0J0bJOS'        ])
			self.plot_DiffMC(self.selections['3J1bJ'    ])
#			self.plot_DiffMC(self.selections['0J0bJOS_loose'  ])
#			self.plot_DiffMC(self.selections['3J1bJOS'    ])
#			self.plot_DiffMC(self.selections['3J0bJOS'    ])
#			self.plot_DiffMC(self.selections['ZElElChMisId'   ])
#			self.plot_DiffMC(self.selections['ZElElChMisId_SS'])
#			self.plot_DiffMC(self.selections['ZMuMu'])
#			self.plot_DiffMC(self.selections['ZElEl'])
#			self.plot_DiffMC(self.selections['ZMuMu_0J'])
#			self.plot_DiffMC(self.selections['ZElEl_0J'])
#			self.plot_DiffMC(self.selections['ZMuMu_0J_noMETCut'])
#			self.plot_DiffMC(self.selections['ZElEl_0J_noMETCut'])


		if DiffPred :
			# produce results tree
			self.skim_tree('Normal')  # makes sure the skimmed SigEvents tree exists
			restree_path = {}
			restree_path['1J0bJ'] = self.path + 'SSDLResults_1J0bJ.root'
			restree_path['2J0bJ'] = self.path + 'SSDLResults_2J0bJ.root'
			restree_path['3J1bJ'] = self.path + 'SSDLResults_3J1bJ.root'
			if '1J0bJ' in restree_path and not os.path.exists(restree_path['1J0bJ']) :
				self.skim_tree('Normal', 1, '1Jskim')  # makes sure the skimmed SigEvents tree exists
				self.make_predictions(self.selections['1J0bJ'], self.path + 'SSDLYields_1Jskim_Normal.root', True)
				if not os.path.exists(restree_path['2J0bJ']) :
					copytree.copytree(restree_path['1J0bJ'], restree_path['2J0bJ'], 'Results', 'NJ > 1')
			if not os.path.exists(restree_path['2J0bJ']) :
				self.make_predictions(self.selections['2J0bJ'], self.path + 'SSDLYields_skim_Normal.root', True)
			if not os.path.exists(restree_path['3J1bJ']) :
				copytree.copytree(restree_path['2J0bJ'], restree_path['3J1bJ'], 'Results', 'NJ > 2 && NbJmed > 0')

			# selections
			sels = {}
#			sels['1J0bJ'    ] = self.selections['1J0bJ'    ]
			sels['1JnobJ_ee'] = self.selections['1JnobJ_ee']
			sels['2J0bJ'    ] = self.selections['2J0bJ'    ]
			sels['2JnobJ_ee'] = self.selections['2JnobJ_ee']
#			sels['3J1bJ'    ] = self.selections['3J1bJ'    ]
#			sels['3J1bJ_ee' ] = self.selections['3J1bJ_ee' ]
			sels['final'    ] = self.selections['final'    ]
#			sels['final++'  ] = self.selections['final++'  ]
#			sels['final--'  ] = self.selections['final--'  ]

			# produce differential predictions
			for name, sel in sels.iteritems() :
				res_path = restree_path['3J1bJ']
				if sel.name.startswith('2J') : res_path = restree_path['2J0bJ']
				if sel.name.startswith('1J') : res_path = restree_path['1J0bJ']
				self.plot_predictions(res_path, sel)


		if IntPred :
			results = self.make_IntPredictions(self.selections['final'], self.path + 'IntPredictions/')
			results_presel = self.make_IntPredictions(self.selections['3J1bJ'], self.path + 'IntPredictions/3J1bJ/', '3J1bJ', noSyst = False)

			# make table of observation and predictions
			tables.make_ObsPredTable(self.path, results['Normal'])

			# make overview table of systematic studies
			tables.make_SystTable(self.path, results, 'al', 'al', self.systematics)

			print '++:   Observed: %3d   Predicted: %5.1f +/- %4.1f' % (
					results['Normal']['++']['al'].obs,
					results['Normal']['++']['al'].tot + results['Normal']['++']['al'].ttw,
					math.sqrt(results['Normal']['++']['al'].tot_err*results['Normal']['++']['al'].tot_err + results['Normal']['++']['al'].ttw_err*results['Normal']['++']['al'].ttw_err))
			print '--:   Observed: %3d   Predicted: %5.1f +/- %4.1f' % (
					results['Normal']['--']['al'].obs,
					results['Normal']['--']['al'].tot + results['Normal']['--']['al'].ttw,
					math.sqrt(results['Normal']['--']['al'].tot_err*results['Normal']['--']['al'].tot_err + results['Normal']['--']['al'].ttw_err*results['Normal']['--']['al'].ttw_err))

			systematics = {}
			systematics['rare'] = self.RareESyst
			tables.make_YieldsTable(self.path, results['Normal']['al'], systematics)
			tables.make_YieldsTable(self.path, results_presel['Normal']['al'], systematics, '3J1bJ')

		if IntMC :
			# selections
			sels = {}
			sels['1J0bJ'    ] = self.selections['1J0bJ'    ]
			sels['2J0bJ'    ] = self.selections['2J0bJ'    ]
#			sels['2JnobJ_ee'] = self.selections['2JnobJ_ee']
			sels['3J1bJ'    ] = self.selections['3J1bJ'    ]
#			sels['3J1bJ_ee' ] = self.selections['3J1bJ_ee' ]
#			sels['final'    ] = self.selections['final'    ]
#			sels['final++'  ] = self.selections['final++'  ]
#			sels['final--'  ] = self.selections['final--'  ]
#			sels['2J0bJOS'  ] = self.selections['2J0bJOS'  ]
#			sels['3J1bJOS'  ] = self.selections['3J1bJOS'  ]
#			sels['4J2bJOS'  ] = self.selections['4J2bJOS'  ]
#			sels['2J0bJOS_ZVeto'] = self.selections['2J0bJOS_ZVeto']
#			sels['2J0bJ_noZVeto'] = self.selections['2J0bJ_noZVeto']

			# produce differential predictions
			for name, sel in sels.iteritems() :
				results = self.make_IntPredictions(sel, self.path + 'IntPredictions/MC/%s' % name, suffix = '_MC_%s' % name, IntMC = True)
				tables.make_MCYieldsTable(self.path, results['al'], suffix = name)

		if FakeClosure :
			# selections
			sels = {}
#			sels['2J0bJ'    ] = self.selections['2J0bJ'    ]
			sels['3J1bJ'    ] = self.selections['3J1bJ'    ]

			ttbartree_path = {}
			ttbartree_path['2J0bJ'] = '%sSSDLYields_TTJets_%s.root' % (self.path, '2J0bJ')
			ttbartree_path['3J1bJ'] = '%sSSDLYields_TTJets_%s.root' % (self.path, '3J1bJ')

			# produce differential predictions
			if not os.path.exists(ttbartree_path['2J0bJ']) :
				copytree.copytree(self.path + 'SSDLYields.root', ttbartree_path['2J0bJ'], 'SigEvents', 'SName == \"TTJets\" && SystFlag == 0 && NJ >= 2')
				if not os.path.exists(ttbartree_path['3J1bJ']) :
					copytree.copytree(ttbartree_path['2J0bJ'], ttbartree_path['3J1bJ'], 'SigEvents', 'SName == \"TTJets\" && SystFlag == 0 && NJ >= 3 && NbJmed >= 1')

			for name, sel in sels.iteritems() :
				closure_res = self.make_closureTest(ttbartree_path[name], ['TTJets'], sel)
#				tables.make_YieldsTable(self.path + 'closure/', closure_res['al'], suffix = name)
				tables.make_closureTable(self.path, closure_res['al'], prefix = 'Fake', suffix = name)

		if SLPlots :
			sl_path = self.path + 'SLYields.root'
			var_sel_str = []
			var_sel_str.append(['ElPFIso', 'NEls == 1'])
			var_sel_str.append(['ElD0'   , 'NEls == 1'])
			var_sel_str.append(['ElDz'   , 'NEls == 1'])
			var_sel_str.append(['MuPFIso', 'NMus == 1'])
			var_sel_str.append(['MuD0'   , 'NMus == 1'])
			var_sel_str.append(['MuDz'   , 'NMus == 1'])
			for [var, sel_str] in var_sel_str :
				self.plot_SL(sl_path, var, sel_str, True)


	def skim_tree(self, syst = '', minNJ = 2, suffix = 'skim') :
		'''skim SigEvents tree'''

		skimtree_path = self.path + 'SSDLYields_%s.root' % suffix
		if not os.path.exists(skimtree_path) :
			print '[status] creating skimmed tree for di-boson, rare MC (SType 15) and data..'
			print '         Data: same-sign, opposite-sign (em, ee), tight-tight, tight-loose, loose-tight, loose-loose events (only nominal, no systematics), NJ >= 2'
			print '         MC:   only same-sign, tight-tight events for all systematics, NJ >= %d' % minNJ
			copytree.copytree(self.path + 'SSDLYields.root', skimtree_path, 'SigEvents', '(SType < 3 || SType == 15) && (SType < 3 || Flavor < 3) && (SType < 3 || TLCat == 0) && (SType > 2 || SystFlag == 0) && (Flavor != 3) && (NJ >= %d)' % minNJ)

		if syst != '' :
			systtree_path = self.path + 'SSDLYields_' + suffix + '_' + syst + '.root'

			if not os.path.exists(systtree_path) :
				print '[status] creating skimmed tree file for %s systematic..' % (syst)
				if syst != 'Normal' : minNJ = 3
				copytree.copytree(skimtree_path, systtree_path, 'SigEvents', 'SystFlag == %d && NJ >= %d' % (self.systflags[syst], minNJ))


	def readDatacard(self, cardfile, verbose = 0, event_count = True) :
		'''
		reads datacard file

		line format:

		data:
		 0: sample name
		 1: input file
		 2: data mc type
		 3: channel (0: DoubleMu, 1: DoubleEle, 2: MuEG, 5: SingleMu)

		mc:
		 0: sample name
		 1: input file
		 2: data mc type
		 3: channel (-1)
		 4: cross section
		'''

		print '[status] reading datacard %s' % (cardfile)
		samples = {}
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
				samples[name] = sample.sample(name = name, datamc = datamc, channel = channel, xsec = xsec, ngen = -1, inputfile = inputfile)
				if verbose > 0 : print samples[name]
		if event_count : self.read_ngen(samples)
		return samples


	def read_ngen(self, samples) :
		for name, s in samples.iteritems() :
			if s.datamc is 0 : continue
			s.ngen = int(self.ssdlfile.Get(s.name+'/'+s.name+'_EventCount').GetEntries())


	def get_samples(self, channel) :
		return sample.sample.get_samples(channel, self.samples)


	def read_selections(self, path) :
		'''reads selections from config file'''

		print '[status] reading selections from %s..' % path
		if not os.path.isfile(path) :
			print '[ERROR] %s does not exist!' % path
			sys.exit(1)
		sel_file = ConfigParser.ConfigParser()
		sel_file.optionxform = str # case sensitive options
		sel_file.read(path)

		sels = {}
		for sel_name in sel_file.sections() :
			sel = selection.selection(name = sel_name)
			for cut in sel_file.options(sel_name) :
				if not hasattr(sel, cut) :
					print '[ERROR] selection instance has no attribute %s! Check your selections config file!' % cut
					sys.exit(1)
				setattr(sel, cut, eval(sel_file.get(sel_name, cut)))
			sels[sel_name] = sel

		return sels


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
		if opt is 0 :
			if charge ==  0 : return 'al'
			if charge == +1 : return '++'
			if charge == -1 : return '--'
			return ''

		elif opt is 1 :
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
		self.ResTree_BECat     = np.zeros(1, dtype=int  ); self.results_tree.Branch('BECat',       self.ResTree_BECat      , 'BECat/I'    );
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


	def get_ChMisID_SF(self) :
		print '[statur] getting charge mis-ID scale factor..'

		self.skim_tree('Normal')  # makes sure the skimmed SigEvents tree exists
		self.chmid_sf = 1.
		sel = copy.deepcopy(self.selections['2JnobJ_ee'])
		sel.mll = 80.
		sel.maxMll = 100.
		res = self.make_predictions(sel, self.path + 'SSDLYields_skim_Normal.root')

		ch_str = 'al'
		chan   = 'al'
		num = res[ch_str][chan].obs - (res[ch_str][chan].tot_exp - res[ch_str][chan].cmid)
		den = res[ch_str][chan].cmid
#		num_err = math.sqrt(res[ch_str][chan].tot_exp_staterr**2 - res[ch_str][chan].cmid_staterr**2)
		num_err = math.sqrt(res[ch_str][chan].tot_exp_err**2 - res[ch_str][chan].cmid_err**2) # take stat+syst error of all other estimations into account
		den_err = res[ch_str][chan].cmid_staterr

		(ratio, ratio_err) = helper.ratio_withError(num, num_err, den, den_err)
		print '         charge mis-ID SF: %4s +- %4s' % helper.get_roundedNumber(ratio, ratio_err)

		return (ratio, ratio_err)


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


	def make_IntPredictions(self, sel, output_path, suffix = '', blind = False, noSyst = False, IntMC = False) :
		'''
		makes predictions for all systematics with a given selection and returns a nested dictionary of result objects
		results[SYST][CHARGE][FLAVOR]
		'''

		if not output_path.endswith('/') :
			output_path += '/'
		results = {}
		resultspath = '%sresults%s.pkl' % (output_path, suffix)

		if os.path.exists(resultspath) :
			print '[status] loading results of predictions from %s..' % (resultspath)
			results = helper.load_object(resultspath)

		elif IntMC :
			results = self.make_IntMC(sel)
			helper.save_object(results, resultspath)

		else :
			for syst in self.systflags :
				if noSyst and syst != 'Normal' : continue
				print '[status] making predictions for %s systematic' % (syst)
				self.skim_tree(syst)  # makes sure the skimmed SigEvents tree exists
				systpath = self.path + 'SSDLYields_skim_' + syst + '.root'
				syst_sel = copy.deepcopy(sel)
				syst_sel.name += '_%s' % syst
				syst_sel.systflag = self.systflags[syst]
				results[syst] = self.make_predictions(syst_sel, systpath, False)
			helper.save_object(results, resultspath)

		if IntMC : return results

		# make datacards for each charge-flavor channel
		datacards_6channels = []
		for charge in results['Normal'] :
			datacards_3channels = []
			for chan in results['Normal'][charge] :
				datacard = self.make_datacard(results, chan, charge, suffix, '%sdatacards%s' % (output_path, suffix), blind)
				ch_str   = results['Normal'][charge][chan].chan_str
				if charge != 'al' and chan != 'al' :
					datacards_6channels.append('%s=%s' % (ch_str, datacard))
				if chan != 'al' :
					datacards_3channels.append('%s=%s' % (ch_str, datacard))
			self.combine_datacards(datacards_3channels, '%sdatacards%s/datacard_ssdl_ttW_3channels_%s%s.txt' % (output_path, suffix, charge, suffix))
		self.combine_datacards(datacards_6channels, '%sdatacards%s/datacard_ssdl_ttW_6channels%s.txt' % (output_path, suffix, suffix))

		return results


	def make_predictions(self, sel, treepath, write_ResTree = False) :
		'''
		oberservation and prediction for different selections

		writes Results tree or returns nested dictionary with result objects
		res[CHARGE][FLAVOR]
		'''

		print '[status] getting observations and predictions with selection %s..' % (sel.name)
		print sel
		if sel.charge != 0 : print '[WARNING] Avoid running make_predictions with a charge selection. Please apply charge selection later.'

		print '[status] open SigEvents tree from %s' % (treepath)
		sigfile = ROOT.TFile.Open(treepath, 'READ')
		sigtree = sigfile.Get('SigEvents')
		print '[status] loaded %s tree with %d events' % (sigtree.GetName(), sigtree.GetEntries())

		##################
		# INIT VARIABLES #
		##################

		FR = fakeratios.FakeRatios()

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
				if charge != 0 or sel.charge != 0 : chargeFactor = 0.5

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
				res[ch_str]['al'].cmid_staterr = math.sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_em_chmid_e1*nt2_em_chmid_e1)
				res[ch_str]['em'].cmid_staterr = nt2_em_chmid_e1
				res[ch_str]['ee'].cmid_staterr = nt2_ee_chmid_e1

		else :
			print '[status] skipped charge mis-ID prediction for systematics'

		################
		# Results tree #
		################

		if write_ResTree :
			print '[status] preparing Results tree to store observations and predictions..'
			self.results_file = ROOT.TFile(self.path + 'SSDLResults_%s.root' % sel.name, 'RECREATE')
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
				self.ResTree_BECat    [0] = event.BECat
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
						res[ch_str]['al'].npp += npp;
						res[ch_str]['al'].npf += npf;
						res[ch_str]['al'].nfp += nfp;
						res[ch_str]['al'].nff += nff;
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

					# apply factor 0.5 as weight and fill ChMisID twice: once for charge = +1 and once for charge = -1
					# in this way it's easier to apply a selection on the charge later
					chargeFactor = 0.5

					# EM OS
					if event.Flavor is 4 :
						if event.BECat is 0 : self.ResTree_Weight[0] = chargeFactor * fbb
						if event.BECat is 1 : self.ResTree_Weight[0] = chargeFactor * fee

					# EE OS
					if event.Flavor is 5 :
						if event.BECat is 0                     : self.ResTree_Weight[0] = chargeFactor * 2*fbb
						if event.BECat is 1 or event.BECat is 2 : self.ResTree_Weight[0] = chargeFactor * 2*feb
						if event.BECat is 3                     : self.ResTree_Weight[0] = chargeFactor * 2*fee

					self.ResTree_Charge[0] = 1
					self.ResTree_CFChan[0] = self.ResTree_Flavor[0]
					self.results_tree.Fill()

					self.ResTree_Charge[0] = -1
					self.ResTree_CFChan[0] = self.ResTree_Flavor[0] + 3
					self.results_tree.Fill()


#				if event.Charge > 0 :
#					self.ResTree_CFChan   [0] = self.ResTree_Flavor[0]
#				else :
#					self.ResTree_CFChan   [0] = self.ResTree_Flavor[0] + 3



#					self.results_tree.Fill()

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

		#####################
		# Fakes predictions #
		#####################

		print '[status] calculating fake predictions..'

		FR.setNToyMCs(100)
		FR.setAddESyst(self.FakeESyst)

		# ratios with ewk subtraction
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
			res[ch_str]['al'].fake_err = math.sqrt(FR.getTotEStat()  *FR.getTotEStat()   + self.FakeESyst2*res[ch_str]['al'].fake*res[ch_str]['al'].fake)
			res[ch_str]['mm'].fake_err = math.sqrt(FR.getMMTotEStat()*FR.getMMTotEStat() + self.FakeESyst2*res[ch_str]['mm'].fake*res[ch_str]['mm'].fake)
			res[ch_str]['em'].fake_err = math.sqrt(FR.getEMTotEStat()*FR.getEMTotEStat() + self.FakeESyst2*res[ch_str]['em'].fake*res[ch_str]['em'].fake)
			res[ch_str]['ee'].fake_err = math.sqrt(FR.getEETotEStat()*FR.getEETotEStat() + self.FakeESyst2*res[ch_str]['ee'].fake*res[ch_str]['ee'].fake)
			res[ch_str]['al'].fake_staterr = FR.getTotEStat()
			res[ch_str]['mm'].fake_staterr = FR.getMMTotEStat()
			res[ch_str]['em'].fake_staterr = FR.getEMTotEStat()
			res[ch_str]['ee'].fake_staterr = FR.getEETotEStat()

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
				res[ch_str][chan].wz           = wz
				res[ch_str][chan].wz_err       = math.sqrt(wz_staterr2 + self.WZESyst2 * wz * wz)
				res[ch_str][chan].wz_staterr   = math.sqrt(wz_staterr2)

				# store ttW mc yields
				res[ch_str][chan].ttw          = ttw
				res[ch_str][chan].ttw_err      = math.sqrt(ttw_staterr2 + self.TTWESyst2 * ttw * ttw)
				res[ch_str][chan].ttw_staterr  = math.sqrt(ttw_staterr2)

				# store ttZ mc yields
				res[ch_str][chan].ttz          = ttz
				res[ch_str][chan].ttz_err      = math.sqrt(ttz_staterr2 + self.TTZESyst2 * ttz * ttz)
				res[ch_str][chan].ttz_staterr  = math.sqrt(ttz_staterr2)

				# store rare mc yields
				res[ch_str][chan].rare         = rare
				res[ch_str][chan].rare_err     = math.sqrt(rare_staterr2 + self.RareESyst2 * rare * rare)
				res[ch_str][chan].rare_staterr = math.sqrt(rare_staterr2)

		# return results only if they are not written to the Results tree
		if not write_ResTree :
			return res


	def make_datacard(self, results, chan, charge, suffix = '', output_dir = '', blind = False) :
		'''
		takes a nested dictionary of result objects as input:
		results[SYSTFLAG][CHARGE][FLAVOR]

		writes datacard in the combine tool format
		'''

		datacard_name = 'datacard_ssdl_ttW_%s%s.txt' % (results['Normal'][charge][chan].chan_str, suffix)
		if output_dir == '' : datacard_path = self.path + 'datacards/'
		else                : datacard_path = '%s/' % (output_dir)
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
			if blind : obs = results['Normal'][charge][chan].ttw + results['Normal'][charge][chan].tot
			else     : obs = results['Normal'][charge][chan].obs
			file.write('observation\t%5.3f\n' % obs)
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
				self.systematics['lumi'],
				self.systematics['lumi'],
				self.systematics['lumi'],
				self.systematics['lumi']))

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
				self.systematics['ltrig'],
				self.systematics['ltrig'],
				self.systematics['ltrig'],
				self.systematics['ltrig']))

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
				self.systematics['scale_up'],
				self.systematics['scale_dn'],
				self.systematics['scale_up'],
				self.systematics['scale_dn'],
				self.systematics['scale_up'],
				self.systematics['scale_dn'],
				self.systematics['scale_up'],
				self.systematics['scale_dn']))

			file.write('tmass    lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t-\t\t-\n' % (
				self.systematics['tmass_up'],
				self.systematics['tmass_dn'],
				self.systematics['tmass_up'],
				self.systematics['tmass_dn']))

			file.write('gen      lnN\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-\n' % (self.systematics['gen']))

			file.write('pdf      lnN\t%5.3f\t\t-\t\t-\t\t-\t\t%5.3f\t\t-\n' % (self.systematics['pdf'], self.systematics['wz_pdf']))

		return datacard_path + datacard_name


	def combine_datacards(self, cards, target_path) :
		'''takes a list of datacards and adds them to a combined datacard'''
		print '[status] combining datacards:\n           %s' % ('\n           '.join(cards))
		print '         target:\n           %s' % target_path
		combineCards = 'combineCards.py %s > %s' % (' '.join(cards), target_path)
		os.system(combineCards)


	def plot_predictions(self, treepath, sel) :
		vars = []
		vars.append('HT'    )
		vars.append('MET'   )
		vars.append('NJ'    )
		vars.append('NbJmed')
		vars.append('pT1'   )
		vars.append('pT2'   )
		vars.append('Mll'   )
		vars.append('Mll_noChMisIDSF')
		vars.append('NVrtx' )
		vars.append('minMT' )
		vars.append('M3'    )
		vars.append('Int'   )
		vars.append('CFChan')  # ChMisID prediction is not totally correct and it's statistical uncertainty a bit too large. Since not all OS events are considered and diveded by 2, but only the ones in the according charge channel.
		vars.append('CFChan_TotalBin')  # ChMisID prediction is not totally correct and it's statistical uncertainty a bit too large. Since not all OS events are considered and diveded by 2, but only the ones in the according charge channel.
		vars.append('Charge')  # ChMisID prediction is not totally correct and it's statistical uncertainty a bit too large. Since not all OS events are considered and diveded by 2, but only the ones in the according charge channel.

		if sel.name == '1JnobJ' or sel.name == '2JnobJ' :
			vars = ['Mll', 'Mll_noChMisIDSF']

		print '[status] open Results tree from %s' % (treepath)
		resfile = ROOT.TFile.Open(treepath, 'READ')
		restree = resfile.Get('Results')

		print '[status] getting plots from %s tree for' % (restree.GetName()), ', '.join(vars)

		for var in vars :
			add_total_bin = False
			if var.endswith('_TotalBin') :
				var = var.split('_')[0]
				add_total_bin = True
			noChMisIDSF = False
			if var.endswith('_noChMisIDSF') :
				var = var.split('_')[0]
				noChMisIDSF = True
			histo_bins = config.get_histoBins(var, sel)
			self.plot_ObsPred(restree, sel, var, histo_bins, add_total_bin, noChMisIDSF)

		resfile.Close()


	def plot_ObsPred(self, tree, sel, var, settings, add_total_bin = False, noChMisIDSF = False) :

		nbins = settings['nbins']
		min   = settings['min'  ]
		max   = settings['max'  ]

		if add_total_bin :
			min = min - (max-min)/nbins
			nbins += 1

		histos = {}
		h_tmp  = {}

		chargeFactor = 1.
		noChargeSelFactor = 0.5
		if sel.charge != 0 or var == 'Charge' or var == 'CFChan' :
			chargeFactor = 0.5
			noChargeSelFactor = 1.

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
		h_wz_name    = 'h_wz_npass'    + var + '_' + sel.name; h_tmp['wz_npass' ] = ROOT.TH1D(h_wz_name   , h_wz_name   , nbins, min, max)
		h_ttz_name   = 'h_ttz_npass'   + var + '_' + sel.name; h_tmp['ttz_npass'] = ROOT.TH1D(h_ttz_name  , h_ttz_name  , nbins, min, max)
		h_ttw_name   = 'h_ttw_npass'   + var + '_' + sel.name; h_tmp['ttw_npass'] = ROOT.TH1D(h_ttw_name  , h_ttw_name  , nbins, min, max)
		tree.Draw(var_str+'>>'+h_wz_name   , '(ObsPred == 3 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_ttz_name  , '(ObsPred == 5 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_ttw_name  , '(ObsPred == 4 && %s)' % sel.get_selectionString(ResTree = True), 'goff')

		# setup histograms for observation, fakes, chmid, wz, ttw, ttz and rare
		h_obs_name   = 'h_obs_'   + var + '_' + sel.name; histos['obs'  ] = ROOT.TH1D(h_obs_name  , h_obs_name  , nbins, min, max)
		h_fake_name  = 'h_fake_'  + var + '_' + sel.name; histos['fake' ] = ROOT.TH1D(h_fake_name , h_fake_name , nbins, min, max)
		h_chmid_name = 'h_chmid_' + var + '_' + sel.name; histos['chmid'] = ROOT.TH1D(h_chmid_name, h_chmid_name, nbins, min, max)
		h_rare_name  = 'h_rare_'  + var + '_' + sel.name; histos['rare' ] = ROOT.TH1D(h_rare_name , h_rare_name , nbins, min, max)
		h_wz_name    = 'h_wz_'    + var + '_' + sel.name; histos['wz'   ] = ROOT.TH1D(h_wz_name   , h_wz_name   , nbins, min, max)
		h_ttz_name   = 'h_ttz_'   + var + '_' + sel.name; histos['ttz'  ] = ROOT.TH1D(h_ttz_name  , h_ttz_name  , nbins, min, max)
		h_ttw_name   = 'h_ttw_'   + var + '_' + sel.name; histos['ttw'  ] = ROOT.TH1D(h_ttw_name  , h_ttw_name  , nbins, min, max)
		h_bgtot_name = 'h_bgtot_' + var + '_' + sel.name; histos['bgtot'] = ROOT.TH1D(h_bgtot_name, h_bgtot_name, nbins, min, max)
		h_pred_name  = 'h_pred_'  + var + '_' + sel.name; histos['pred' ] = ROOT.TH1D(h_pred_name , h_pred_name , nbins, min, max)
		h_mm_nt2_name  = 'h_mm_nt2_'  + var + '_' + sel.name; h_tmp['mm_nt2_npass' ] = ROOT.TH1D(h_mm_nt2_name , h_mm_nt2_name , nbins, min, max)
		h_mm_nt10_name = 'h_mm_nt10_' + var + '_' + sel.name; h_tmp['mm_nt10_npass'] = ROOT.TH1D(h_mm_nt10_name, h_mm_nt10_name, nbins, min, max)
		h_mm_nt01_name = 'h_mm_nt01_' + var + '_' + sel.name; h_tmp['mm_nt01_npass'] = ROOT.TH1D(h_mm_nt01_name, h_mm_nt01_name, nbins, min, max)
		h_mm_nt0_name  = 'h_mm_nt0_'  + var + '_' + sel.name; h_tmp['mm_nt0_npass' ] = ROOT.TH1D(h_mm_nt0_name , h_mm_nt0_name , nbins, min, max)
		h_em_nt2_name  = 'h_em_nt2_'  + var + '_' + sel.name; h_tmp['em_nt2_npass' ] = ROOT.TH1D(h_em_nt2_name , h_em_nt2_name , nbins, min, max)
		h_em_nt10_name = 'h_em_nt10_' + var + '_' + sel.name; h_tmp['em_nt10_npass'] = ROOT.TH1D(h_em_nt10_name, h_em_nt10_name, nbins, min, max)
		h_em_nt01_name = 'h_em_nt01_' + var + '_' + sel.name; h_tmp['em_nt01_npass'] = ROOT.TH1D(h_em_nt01_name, h_em_nt01_name, nbins, min, max)
		h_em_nt0_name  = 'h_em_nt0_'  + var + '_' + sel.name; h_tmp['em_nt0_npass' ] = ROOT.TH1D(h_em_nt0_name , h_em_nt0_name , nbins, min, max)
		h_ee_nt2_name  = 'h_ee_nt2_'  + var + '_' + sel.name; h_tmp['ee_nt2_npass' ] = ROOT.TH1D(h_ee_nt2_name , h_ee_nt2_name , nbins, min, max)
		h_ee_nt10_name = 'h_ee_nt10_' + var + '_' + sel.name; h_tmp['ee_nt10_npass'] = ROOT.TH1D(h_ee_nt10_name, h_ee_nt10_name, nbins, min, max)
		h_ee_nt01_name = 'h_ee_nt01_' + var + '_' + sel.name; h_tmp['ee_nt01_npass'] = ROOT.TH1D(h_ee_nt01_name, h_ee_nt01_name, nbins, min, max)
		h_ee_nt0_name  = 'h_ee_nt0_'  + var + '_' + sel.name; h_tmp['ee_nt0_npass' ] = ROOT.TH1D(h_ee_nt0_name , h_ee_nt0_name , nbins, min, max)
		h_nt2_em_BB_os_name = 'h_nt2_em_BB_os_' + var + '_' + sel.name; h_tmp['nt2_em_BB_os'] = ROOT.TH1D(h_nt2_em_BB_os_name, h_nt2_em_BB_os_name, nbins, min, max)
		h_nt2_em_EE_os_name = 'h_nt2_em_EE_os_' + var + '_' + sel.name; h_tmp['nt2_em_EE_os'] = ROOT.TH1D(h_nt2_em_EE_os_name, h_nt2_em_EE_os_name, nbins, min, max)
		h_nt2_ee_BB_os_name = 'h_nt2_ee_BB_os_' + var + '_' + sel.name; h_tmp['nt2_ee_BB_os'] = ROOT.TH1D(h_nt2_ee_BB_os_name, h_nt2_ee_BB_os_name, nbins, min, max)
		h_nt2_ee_EB_os_name = 'h_nt2_ee_EB_os_' + var + '_' + sel.name; h_tmp['nt2_ee_EB_os'] = ROOT.TH1D(h_nt2_ee_EB_os_name, h_nt2_ee_EB_os_name, nbins, min, max)
		h_nt2_ee_EE_os_name = 'h_nt2_ee_EE_os_' + var + '_' + sel.name; h_tmp['nt2_ee_EE_os'] = ROOT.TH1D(h_nt2_ee_EE_os_name, h_nt2_ee_EE_os_name, nbins, min, max)

		# getting histograms from results tree
		tree.Draw(var_str+'>>'+h_obs_name  ,    'Weight*(ObsPred == 0 && TLCat == 0 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_fake_name ,    'Weight*(ObsPred == 1 && %s)'               % sel.get_selectionString(ResTree = True), 'goff')
		tree.Draw(var_str+'>>'+h_chmid_name,    'Weight*(ObsPred == 2 && TLCat == 0 && %s)' % sel.get_selectionString(ResTree = True), 'goff')
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
		tree.Draw(var_str+'>>'+h_nt2_em_BB_os_name, '%f*(ObsPred == 2 && TLCat == 0 && BECat == 0 && Flavor == 1 && %s)' % (noChargeSelFactor, sel.get_selectionString(ResTree = True)), 'goff')
		tree.Draw(var_str+'>>'+h_nt2_em_EE_os_name, '%f*(ObsPred == 2 && TLCat == 0 && BECat == 1 && Flavor == 1 && %s)' % (noChargeSelFactor, sel.get_selectionString(ResTree = True)), 'goff')
		tree.Draw(var_str+'>>'+h_nt2_ee_BB_os_name, '%f*(ObsPred == 2 && TLCat == 0 && BECat == 0 && Flavor == 2 && %s)' % (noChargeSelFactor, sel.get_selectionString(ResTree = True)), 'goff')
		tree.Draw(var_str+'>>'+h_nt2_ee_EB_os_name, '%f*(ObsPred == 2 && TLCat == 0 && (BECat == 1 || BECat == 2) && Flavor == 2 && %s)' % (noChargeSelFactor, sel.get_selectionString(ResTree = True)), 'goff')
		tree.Draw(var_str+'>>'+h_nt2_ee_EE_os_name, '%f*(ObsPred == 2 && TLCat == 0 && BECat == 3 && Flavor == 2 && %s)' % (noChargeSelFactor, sel.get_selectionString(ResTree = True)), 'goff')

		if add_total_bin :
			for key, histo in dict(histos.items() + h_tmp.items() + h_rares_npass.items()).iteritems() :
				sum = histo.Integral()
				histo.SetBinContent(1, sum)

		if noChMisIDSF :
			histos['chmid'].Scale(1./self.chmid_sf)
			chmid_sf = 1.
		else :
			chmid_sf = self.chmid_sf

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

		FR = fakeratios.FakeRatios()
		FR.setNToyMCs(100)
		FR.setAddESyst(self.FakeESyst)

		FR.setMFRatio(self.fpr.MufRatio, self.fpr.MufRatioE) # set error to pure statistical of ratio
		FR.setEFRatio(self.fpr.ElfRatio, self.fpr.ElfRatioE)
		FR.setMPRatio(self.fpr.MupRatio, self.fpr.MupRatioE)
		FR.setEPRatio(self.fpr.ElpRatio, self.fpr.ElpRatioE)

		(fbb, fbbE) = self.calculateChMisIdProb(self.get_samples('DoubleEle'), 'BB', chmid_sf)
		(feb, febE) = self.calculateChMisIdProb(self.get_samples('DoubleEle'), 'EB', chmid_sf)
		(fee, feeE) = self.calculateChMisIdProb(self.get_samples('DoubleEle'), 'EE', chmid_sf)

		for bin in range(1, histos['bgtot'].GetNbinsX()+1) :
			err2 = 0.

			# fakes
			FR.setMMNtl(h_tmp['mm_nt2_npass'].GetBinContent(bin), h_tmp['mm_nt10_npass'].GetBinContent(bin), h_tmp['mm_nt01_npass'].GetBinContent(bin), h_tmp['mm_nt0_npass'].GetBinContent(bin))
			FR.setEMNtl(h_tmp['em_nt2_npass'].GetBinContent(bin), h_tmp['em_nt10_npass'].GetBinContent(bin), h_tmp['em_nt01_npass'].GetBinContent(bin), h_tmp['em_nt0_npass'].GetBinContent(bin))
			FR.setEENtl(h_tmp['ee_nt2_npass'].GetBinContent(bin), h_tmp['ee_nt10_npass'].GetBinContent(bin), h_tmp['ee_nt01_npass'].GetBinContent(bin), h_tmp['ee_nt0_npass'].GetBinContent(bin))

			fake_nPass = histos['fake'].GetBinContent(bin)
			fake_syst2 = self.FakeESyst2 * fake_nPass*fake_nPass
			fake_stat2 = FR.getTotEStat() * FR.getTotEStat()
			err2 += fake_syst2 + fake_stat2
#			print 'bin %d: %f +/- %f (syst) +/- %f (stat)' % (bin, fake_nPass, math.sqrt(fake_syst2), math.sqrt(fake_stat2))

			# charge mis ID
			nt2_em_BB_os = h_tmp['nt2_em_BB_os'].GetBinContent(bin)
			nt2_em_EE_os = h_tmp['nt2_em_EE_os'].GetBinContent(bin)
			nt2_ee_BB_os = h_tmp['nt2_ee_BB_os'].GetBinContent(bin)
			nt2_ee_EB_os = h_tmp['nt2_ee_EB_os'].GetBinContent(bin)
			nt2_ee_EE_os = h_tmp['nt2_ee_EE_os'].GetBinContent(bin)

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
			wz_nPass = h_tmp['wz_npass'].GetBinContent(bin)
			wz_stat2 = scale*scale * self.samples[wz_str].getError2(wz_nPass)
			wz_yield = histos['wz'].GetBinContent(bin)
			wz_syst2 = self.WZESyst2 * wz_yield*wz_yield
			err2 += wz_syst2 + wz_stat2

			# ttz
			scale = self.lumi / self.samples[ttz_str].getLumi()
			ttz_nPass = h_tmp['ttz_npass'].GetBinContent(bin)
			ttz_stat2 = scale*scale * self.samples[ttz_str].getError2(ttz_nPass)
			ttz_yield = histos['ttz'].GetBinContent(bin)
			ttz_syst2 = self.TTZESyst2 * ttz_yield*ttz_yield
			err2 += ttz_syst2 + ttz_stat2

#			# ttz
#			scale = self.lumi / self.samples[ttw_str].getLumi()
#			ttw_nPass = h_tmp['ttw_npass'].GetBinContent(bin)
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

		histofile = ROOT.TFile.Open(self.path + 'SSDLHistos.root', 'UPDATE')
		histofile.cd()
		if not histofile.cd('%s' % sel.name) :
			histofile.mkdir('%s' % sel.name)
			histofile.cd   ('%s' % sel.name)
#		histofile.mkdir('%s' % var)
#		histofile.cd   ('%s' % var)
		for process, histo in histos.iteritems() :
			histo.Write()
		histofile.Close()

		if add_total_bin : var += '_TotalBin'
		cms_label = 0
		if   sel.flavor == -2 : prefix = '_OS'
		elif sel.flavor == -1 : prefix = '_ALL'
		elif sel.flavor ==  0 : prefix = '_MM'
		elif sel.flavor ==  1 : prefix = '_EM'
		elif sel.flavor ==  2 : prefix = '_EE'
		else                  : prefix = ''
		var_tabstr = var
		if var == 'Mll' :
			if   sel.flavor == 0 or sel.flavor == 3 : var_tabstr = 'Mmm'
			elif sel.flavor == 1 or sel.flavor == 4 : var_tabstr = 'Mem'
			elif sel.flavor == 2 or sel.flavor == 5 : var_tabstr = 'Mee'
		charge_str = ''
		if sel.charge != 0 :
			if sel.charge > 0 : charge_str = '^{+}'
			if sel.charge < 0 : charge_str = '^{-}'
		if sel.name == 'final' : cms_label = 0
#		cms_label = 2

		# save histo data to table
		pl = ttvStyle.ttvStyle()
		suffix = ''
		if noChMisIDSF : suffix = '_noChMisIDSF'
		path = '%sObsPredPlots/%s/' % (self.path, sel.name)
		helper.mkdir(path)
		helper.save_histo2table(histos = histos, processes = ['obs', 'pred', 'fake', 'rare', 'chmid', 'wz', 'ttz', 'ttw'], path = '%sObsPred%s_%s%s.dat' % (path, prefix, var, suffix), var = var_tabstr, lumi = self.lumi, bin_width = pl.label_binWidth(var), asymmErr = ['obs'])

		# save histo
		for TeX_switch in [True, False] :
			pl = ttvplot.ttvplot(self.path + 'ObsPredPlots/%s/'%sel.name, '2L', lumi = self.lumi, cms_label = cms_label, asymmErr = True, TeX_switch = TeX_switch)
			pl.save_plot(histos, var, prefix = prefix, suffix = suffix, charge_str = charge_str)
			# preliminary plots with suffix
			pl.cms_label = 2
#			pl.save_plot(histos, var, prefix = prefix, suffix = 'preliminary', charge_str = charge_str)
#		raw_input('ok? ')
#		raw_input('ok? ')


	def make_IntMC(self, sel) :
		'''
		loops over SigEvents tree and stores all MC yields for a given selection in a nested dictionary of result objects
		res[CHARGE][CHANNEL]
		'''

		if sel.flavor < -1 or sel.flavor > 2 :
			print '[WARNING] Yields split by charge doe not make sense for an opposite-sign selection!'
		skimtree_path = '%sSSDLYields_%s.root' % (self.path, sel.name)
		if not os.path.exists(skimtree_path) :
			copytree.copytree('%sSSDLYields.root' % self.path, skimtree_path, 'SigEvents', sel.get_selectionString())
		file = ROOT.TFile.Open(skimtree_path, 'READ')
		tree = file.Get('SigEvents')
		last_sample = ''
		mc       = {}
		mc_npass = {}
		res      = {}
		missing_samples = []
		for ch_str, charge in self.charges.iteritems() :
			res[ch_str] = {}
			charge_str = ''
			if charge > 0 : charge_str = '+'
			if charge < 0 : charge_str = '-'
			res[ch_str]['al'] = result.result('al', charge, 'int'+charge_str+charge_str)
			res[ch_str]['mm'] = result.result('mm', charge, 'm'+charge_str+'m'+charge_str)
			res[ch_str]['em'] = result.result('em', charge, 'e'+charge_str+'m'+charge_str)
			res[ch_str]['ee'] = result.result('ee', charge, 'e'+charge_str+'e'+charge_str)
		for s in self.samples :
			mc      [s] = {}
			mc_npass[s] = {}
			for ch_str, charge in self.charges.iteritems() :
				mc      [s][ch_str] = {}
				mc_npass[s][ch_str] = {}
				for chan in res[ch_str] :
					mc      [s][ch_str][chan] = 0.
					mc_npass[s][ch_str][chan] = 0

		print '[status] looping over SigEvents tree from %s with %d events..' % (skimtree_path, tree.GetEntries())
		for event in tree :
			if last_sample != str(event.SName) :
				print '[status] processing %s..' % (event.SName)
				last_sample = str(event.SName)

			if not sel.passes_selection(event = event) : continue

			if str(event.SName) not in self.samples :
				if str(event.SName) not in missing_samples :
					print '[WARNING] %s was not found!' % str(event.SName)
					missing_samples.append(str(event.SName))
				continue

			# GET ALL DATA EVENTS
			if event.SType < 3 :
				for ch_str, charge in self.charges.iteritems() :
					if charge != 0 and event.Charge != charge : continue

					# int
					if event.TLCat is 0 : res[ch_str]['al'].nt2  += 1
					if event.TLCat is 1 : res[ch_str]['al'].nt10 += 1
					if event.TLCat is 2 : res[ch_str]['al'].nt01 += 1
					if event.TLCat is 3 : res[ch_str]['al'].nt0  += 1

					# MM
					if event.Flavor is 0 or event.Flavor is 3 :
						if event.TLCat is 0 : res[ch_str]['mm'].nt2  += 1
						if event.TLCat is 1 : res[ch_str]['mm'].nt10 += 1
						if event.TLCat is 2 : res[ch_str]['mm'].nt01 += 1
						if event.TLCat is 3 : res[ch_str]['mm'].nt0  += 1

					# EM
					if event.Flavor is 1 or event.Flavor is 4 :
						if event.TLCat is 0 : res[ch_str]['em'].nt2  += 1
						if event.TLCat is 1 : res[ch_str]['em'].nt10 += 1
						if event.TLCat is 2 : res[ch_str]['em'].nt01 += 1
						if event.TLCat is 3 : res[ch_str]['em'].nt0  += 1

					# EE
					if event.Flavor is 2 or event.Flavor is 5 :
						if event.TLCat is 0 : res[ch_str]['ee'].nt2  += 1
						if event.TLCat is 1 : res[ch_str]['ee'].nt10 += 1
						if event.TLCat is 2 : res[ch_str]['ee'].nt01 += 1
						if event.TLCat is 3 : res[ch_str]['ee'].nt0  += 1

			# GET MC EVENTS
			if event.SType > 5 and event.TLCat == 0 :
				scale = event.PUWeight * event.HLTSF * self.lumi / self.samples[str(event.SName)].getLumi()

				if str(event.SName) not in mc :
					mc      [str(event.SName)] = {}
					mc_npass[str(event.SName)] = {}
					for ch_str, charge in self.charges.iteritems() :
						mc      [str(event.SName)][ch_str] = {}
						mc_npass[str(event.SName)][ch_str] = {}
						for chan in res[ch_str] :
							mc      [str(event.SName)][ch_str][chan] = 0.
							mc_npass[str(event.SName)][ch_str][chan] = 0

				for ch_str, charge in self.charges.iteritems() :
					if charge != 0 and event.Charge != charge : continue

					mc      [str(event.SName)][ch_str]['al'] += scale
					mc_npass[str(event.SName)][ch_str]['al'] += 1

					if event.Flavor is 0 or event.Flavor is 3 :
						mc      [str(event.SName)][ch_str]['mm'] += scale
						mc_npass[str(event.SName)][ch_str]['mm'] += 1

					if event.Flavor is 1 or event.Flavor is 4 :
						mc      [str(event.SName)][ch_str]['em'] += scale
						mc_npass[str(event.SName)][ch_str]['em'] += 1

					if event.Flavor is 2 or event.Flavor is 5 :
						mc      [str(event.SName)][ch_str]['ee'] += scale
						mc_npass[str(event.SName)][ch_str]['ee'] += 1

		for ch_str in res :
			for chan in res[ch_str] :
#				if ch_str != 'al' or chan != 'al' : continue

				for s in mc :
					scale = self.lumi / self.samples[s].getLumi()
					staterr = scale * self.samples[s].getError(mc_npass[s][ch_str][chan])

					if s in self.get_samples('Top') :
						if s.startswith('Single') :
							#print 'single top: %s' % s
							sname = '\\singletop'
						else :
							#print 'ttbar: %s' % s
							sname = '\\ttbar'

					elif s in self.get_samples('WJets') :
						#print 'wjets: %s' % s
						sname = '\\wjets'

					elif s in ['DYJets', 'DYJets10To50'] :
						#print 'zjets: %s' % s
						sname = '\\zjets'

					elif s.startswith('GJets') :
						#print 'gjets: %s' % s
						sname = '\\gjets'

					elif s == 'WWTo2L2Nu' :
						sname = '\\ww'

					elif s == 'WZTo3LNu' :
						sname = '\\wz'

					elif s == 'ZZTo4L' :
						sname = '\\zz'

					elif s == 'TTbarW' :
						sname = '\\ttw'

					elif s == 'TTbarZ' :
						sname = '\\ttz'

					elif s in self.get_samples('TTH') :
						sname = '\\ttH'

					elif s in self.get_samples('Multiboson') :
						sname = '\\triboson'

					elif s in ['W+W+', 'W-W-'] :
						#print 'W+-W+-: %s' % s
						sname = '\\WpmWpm'

					elif s in ['WGstarTau', 'WGstarMu'] :
						sname = '\\WGstar'

					elif s == 'WbbJets' :
						sname = '\\WbbJets'

					elif s == 'TTbarG' :
						sname = '\\TTbarG'

					elif s == 'TTbarWW' :
						sname = '\\TTbarWW'

					elif s == 'DPSWW' :
						sname = '\\DPSWW'

					elif s == 'TbZ' :
						sname = '\\TbZ'

					else :
						sname = s

					for name in [sname, 'tot', 'bg'] :
						if name == 'bg' and s == 'TTbarW' : continue
						if name not in res[ch_str][chan].mc :
							res[ch_str][chan].mc[name]         = 0.
							res[ch_str][chan].mc_staterr[name] = 0.
							res[ch_str][chan].mc_xsec[name]    = 0.

						res[ch_str][chan].mc[name]         += mc[s][ch_str][chan]
						res[ch_str][chan].mc_staterr[name]  = math.sqrt(res[ch_str][chan].mc_staterr[name]**2 + staterr**2)
						res[ch_str][chan].mc_xsec[name]    += pl.samples[s].xsec

		return res


	def plot_DiffMC(self, sel) :
		skimtree_path = '%sSSDLYields_%s.root' % (self.path, sel.name)
		if not os.path.exists(skimtree_path) :
			copytree.copytree('%sSSDLYields.root' % self.path, skimtree_path, 'SigEvents', sel.get_selectionString())
		file = ROOT.TFile.Open(skimtree_path, 'READ')
		tree = file.Get('SigEvents')
		if sel.minNjets < 1 :
			self.plot_ObsMC(tree, sel, 'NVrtx', config.get_histoBins('NVrtx', sel), weight_str = 'LepSF', plot_shapes = False, pu_weight = True)
			self.plot_ObsMC(tree, sel, 'NVrtx', config.get_histoBins('NVrtx', sel), weight_str = 'LepSF', plot_shapes = False, pu_weight = False)
			return
		vars = []
#		vars.append('Mll'     )
		vars.append('pT2'     )
		vars.append('NbJmed'  )
		vars.append('HT'      )
#		vars.append('NVrtx'   )
		vars.append('PFIso1'  )
#		vars.append('PFIso2'  )
		vars.append('NJ'      )
		vars.append('MET'     )
		for var in vars :
			settings = config.get_histoBins(var, sel)
			if var == 'HT' :
				settings['nbins'] *= 2
			self.plot_ObsMC(tree, sel, var, settings, plot_shapes = ['ttw', 'ttbar'], pu_weight = True)
#			self.plot_ObsMC(tree, sel, var, settings, plot_shapes = ['ttbar', 'qcd'], pu_weight = True)


	def plot_ObsMC(self, tree, sel, var, settings, add_total_bin = False, weight_str = 'HLTSF', pu_weight = True, plot_shapes = False, plot_missing_samples = False, only_shapes = False) :

		path = '%sObsMCPlots/%s/' % (self.path, sel.name)
		nbins = settings['nbins']
		min   = settings['min'  ]
		max   = settings['max'  ]
		if add_total_bin :
			min = min - (max-min)/nbins
			nbins += 1

		histos = {}
		shapes = {}

		if pu_weight : weight_str += '*PUWeight'

		############################
		# SETUP AND GET HISTOGRAMS #
		############################

		if   var == 'pT1'   : var_str = 'TMath::Max(pT1,pT2)'
		elif var == 'pT2'   : var_str = 'TMath::Min(pT1,pT2)'
		elif var == 'minMT' : var_str = 'TMath::Min(MTLep1,MTLep2)'
		elif var == 'Int'   : var_str = 'Flavor'
		else                : var_str = var

		if not only_shapes :

			# processes
			processes = []
			processes.append('top'  )
			processes.append('zjets')
			processes.append('wjets')
			processes.append('rare' )
			processes.append('wz'   )
			processes.append('qcd'  )
			processes.append('ttz'  )
			processes.append('ttw'  )

			# setup histograms for observation, fakes, chmid, wz, ttw, ttz and rare
			for process in processes :
				h_name   = 'h_%s_%s_%s' % (process, var, sel.name)
				histos[process] = self.get_mcHistoFromTree(tree, self.get_samples(process), var_str, h_name, settings, weight_str, sel.get_selectionString())
			h_obs_name   = 'h_obs_'   + var + sel.name; histos['obs'  ] = ROOT.TH1D(h_obs_name  , h_obs_name  , nbins, min, max)
			h_bgtot_name = 'h_bgtot_' + var + sel.name; histos['bgtot'] = ROOT.TH1D(h_bgtot_name, h_bgtot_name, nbins, min, max)
			h_pred_name  = 'h_pred_'  + var + sel.name; histos['pred' ] = ROOT.TH1D(h_pred_name , h_pred_name , nbins, min, max); histos['pred' ].Sumw2()
			h_stack_name = 'h_stack_' + var + sel.name; histos['stack'] = ROOT.THStack(h_stack_name, h_stack_name)

			# missing samples
			missing = []
			for sample in self.samples :
				if any([sample in self.get_samples(i) for i in processes]) : continue
				if sample in self.get_samples('DoubleMu' ) : continue
				if sample in self.get_samples('DoubleEle') : continue
				if sample in self.get_samples('MuEG'     ) : continue
				if sample in self.get_samples('SingleMu' ) : continue
				if sample in ['WJets2', 'WJets1'] : continue
				missing.append(sample)
			print '[info] NOT considered samples: %s' % ', '.join(missing)
			if plot_missing_samples :
				h_miss_name  = 'h_miss_'  + var + sel.name; histos['miss' ] = self.get_mcHistoFromTree(tree, missing                   , var_str, h_miss_name , settings, weight_str, sel.get_selectionString())
				processes.insert(0, 'miss')

			# getting data
			print '[status] getting %s histogram from data..' % var
			tree.Draw(var_str+'>>'+h_obs_name, 'SType < 3 && %s' % sel.get_selectionString(), 'goff')

			# adding mc samples
			for histo in [histos['bgtot'], histos['pred'], histos['stack']] :
				for process in processes :
					if process not in histos : continue
					if 'h_bgtot' in histo.GetName() and process == 'ttw' : continue
					histo.Add(histos[process])

			# save histo data to table
			prefix = ''
			suffix = ''
			if not pu_weight : suffix = '_noPUWeight'
			if var == 'NVrtx' : bin_width = False
			else              : bin_width = True
			helper.mkdir(path)
			helper.save_histo2table(histos, processes+['obs', 'pred'], '%sObsMC%s_%s%s.dat' % (path, prefix, var, suffix), var, self.lumi, bin_width)

			set_maximum = True

			for TeX_switch in [True, False] :
				pl = ttvStyle.ttvStyle(lumi = self.lumi, cms_label = 0, TeX_switch = TeX_switch)
				canvas = pl.get_canvas(var)
				canvas.cd()

				# legend
				leg_entries = []
				leg_entries.append([histos['obs'], pl.get_processName('obs'), 'lp'])
				for process in reversed(processes) :
					if process not in histos : continue
					histo = histos[process]
					leg_entries.append([histo, pl.get_processName(process), 'f'])
					histo.SetFillColor(pl.get_fillColor(process))
				leg_entries.append([histos['pred'], 'MC uncertainty', 'fl'])
				leg = pl.draw_legend(leg_entries)

				if set_maximum :
					pl.get_maximum(histos.values())
					set_maximum = False

				histos['pred'].SetLineWidth(0)
				histos['pred'].SetMarkerSize(0)
				histos['pred'].SetFillColor(12)
				histos['pred'].SetFillStyle(3005)

#				histos['obs'].SetMarkerStyle(20)
#				histos['pred'].Draw()
				histos['stack'].Draw('hist')
				histos['stack'].GetXaxis().SetTitle(pl.get_varName(var))
				histos['stack'].GetYaxis().SetTitle('Events')
				leg.Draw()
				histos['pred'].Draw('0 E2 same')
#				histos['bgtot'].Draw('hist same')
				histos['obs'  ].Draw('PE X0 same')
				pl.draw_cmsLine()

				canvas.Update()
				if TeX_switch : format_str = 'tex'
				else          : format_str = 'pdf'
				canvas.Print('%sObsMC%s_%s%s.%s' % (path, prefix, var, suffix, format_str))
				ROOT.gPad.SetLogy()
				canvas.Update()
				canvas.Print('%sObsMC%s_%s%s_log.%s' % (path, prefix, var, suffix, format_str))

		if plot_shapes != False :
		# histograms for shape plots
			for process in plot_shapes :
				if process in histos :
					shapes[process] = histos[process].Clone()
				else :
					h_name   = 'h_%s_%s_%s' % (process, var, sel.name)
					shapes[process] = self.get_mcHistoFromTree(tree, self.get_samples(process), var_str, h_name, settings, weight_str, sel.get_selectionString())

			for TeX_switch in [True, False] :
				pl = ttvStyle.ttvStyle(lumi = self.lumi, cms_label = 0, TeX_switch = TeX_switch)
				canvas = pl.get_canvas(var)
				canvas.cd()

				pl.cms_label = 1
				ROOT.gPad.SetLogy(0)
				leg_entries = []

				for process in plot_shapes :
					shape = shapes[process]
					shape.Scale(1. / shape.Integral())
					leg_entries.append([shape, pl.get_processName(process), 'f'])
					shape.SetLineColor(pl.get_fillColor(process))
					shape.SetFillColor(pl.get_fillColor(process))
					shape.SetFillStyle(pl.get_fillStyle(process))

				leg = pl.draw_legend(leg_entries)
				pl.get_maximum(shapes.values())

				for index, process in enumerate(plot_shapes) :
					if index == 0 : draw_option = 'hist'
					else          : draw_option = 'hist same'
					shapes[process].Draw(draw_option)

				shapes[plot_shapes[0]].GetXaxis().SetTitle(pl.get_varName(var))
				shapes[plot_shapes[0]].GetYaxis().SetTitle('Normalised to Unity')

				pl.draw_cmsLine()
				leg.Draw()

				shapes_path = '%sShapes/' % path
				prefix = ''
				suffix = ''
				if not pu_weight : suffix = '_noPUWeight'
				helper.mkdir(shapes_path)
				canvas.Update()
				if TeX_switch : format_str = 'tex'
				else          : format_str = 'pdf'
				canvas.Print('%sObsMC%s_%s%s.%s' % (shapes_path, prefix, var, suffix, format_str))
				helper.save_histo2table(histos = shapes, processes = plot_shapes, path = '%sObsMC%s_%s%s.dat' % (shapes_path, prefix, var, suffix), var = var, lumi = self.lumi, bin_width = False, last_bin = True)
				ROOT.gPad.SetLogy()
				canvas.Update()
				canvas.Print('%sObsMC%s_%s%s_log.%s' % (shapes_path, prefix, var, suffix, format_str))


	def plot_SL(self, sl_path, var, sel_str, ttbar_true = False) :
		sl_file = ROOT.TFile.Open(sl_path, 'READ')
		tree = sl_file.Get('SLEvents')
		path = '%sSLPlots/Shapes/' % self.path
		helper.mkdir(path)

		# processes
		processes = ['ttbar', 'qcd']

		# get histograms from tree
		histos = {}
		scale = {}
		settings = config.get_histoBins(var)
		for process in processes :
			h_name = 'h_%s_%s' % (process, var)
			histos[process] = self.get_mcHistoFromTree(tree = tree, samples = self.get_samples(process), var = var, name = h_name, settings = settings, sel_str = sel_str)
			scale[process] = 1./histos[process].Integral()
			histos[process].Scale(scale[process])

		if ttbar_true :
			if 'El' in var : sel_str += ' && IsSignalElectron == 1'
			if 'Mu' in var : sel_str += ' && IsSignalMuon == 1'

			process = 'ttbar_true'
			h_name  = 'h_ttbar_true_%s' % var
			histos[process] = self.get_mcHistoFromTree(tree = tree, samples = self.get_samples('ttbar'), var = var, name = h_name, settings = settings, sel_str = sel_str)
			histos[process].Scale(scale['ttbar'])
			processes.append(process)

		helper.save_histo2table(histos = histos, processes = processes, path = '%sSLPlot_%s.dat' % (path, var), var = var)

		sl_file.Close()


	def get_mcHistoFromTree(self, tree, samples, var, name, settings, weight = '1.', sel_str = '1==1') :
		'''getting histogram for a list of samples'''

		histo = ROOT.TH1D(name, name, settings['nbins'], settings['min'], settings['max'])
		histo.Sumw2()
		if not sel_str.startswith('&&') : sel_str = '&& ' + sel_str
		for i, sample in enumerate(samples) :
			sys.stdout.write('[status] getting %s histogram from %-16s' % (var, sample + '..'))
			sys.stdout.flush()
			h_tmp_name = 'h_tmp_%d_%d' % (i, self.rand.Integer(10000))
			h_tmp = ROOT.TH1D(h_tmp_name, h_tmp_name, settings['nbins'], settings['min'], settings['max'])
			h_tmp.Sumw2()
			tree.Draw(var+'>>'+h_tmp_name, '%s*(SName == \"%s\" %s)' % (weight, sample, sel_str), 'goff')
			histo.Add(h_tmp, self.lumi / self.samples[sample].getLumi())
			print '%10d events: %10.1f * %8.2f = %10.1f' % (h_tmp.GetEntries(), h_tmp.Integral(), self.lumi / self.samples[sample].getLumi(), self.lumi / self.samples[sample].getLumi() * h_tmp.Integral())
		return histo


	def make_closureTest(self, tree_path, samples, sel) :

		print tree_path
		file = ROOT.TFile.Open(tree_path, 'READ')
		tree = file.Get('SigEvents')

		##################
		# INIT VARIABLES #
		##################

		FR = fakeratios.FakeRatios()

		nt2_npass = {}
		nt10_npass = {}
		nt01_npass = {}
		nt0_npass = {}
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

		last_sample = ''

		for event in tree :
			if last_sample != str(event.SName) :
				print '[status] processing %s..' % (event.SName)
				last_sample = str(event.SName)

			if not sel.passes_selection(event = event, ttLeptons = False, OSwoZVeto = True) : continue
			chan = self.get_channelString(int(event.Flavor))

			if str(event.SName) not in samples : continue
			if event.Flavor > 2                : continue

			scale = self.lumi / self.samples[str(event.SName)].getLumi()# * event.PUWeight * event.HLTSF
			#scale = 1.

			# fake, prompt predictions
			if chan is 'ElMu' :
				f1 = self.fpr.get_fRatio('Muon', event.pT1, event.eta1, 1)
				f2 = self.fpr.get_fRatio('Elec', event.pT2, event.eta2, 1)
				p1 = self.fpr.get_pRatio('Muon', event.pT1, 1)
				p2 = self.fpr.get_pRatio('Elec', event.pT2, 1)
			else :
				f1 = self.fpr.get_fRatio(chan, event.pT1, event.eta1, 1)
				f2 = self.fpr.get_fRatio(chan, event.pT2, event.eta2, 1)
				p1 = self.fpr.get_pRatio(chan, event.pT1, 0)
				p2 = self.fpr.get_pRatio(chan, event.pT2, 0)

			npp = scale * FR.getWpp(event.TLCat, f1, f2, p1, p2)
			npf = scale * FR.getWpf(event.TLCat, f1, f2, p1, p2)
			nfp = scale * FR.getWfp(event.TLCat, f1, f2, p1, p2)
			nff = scale * FR.getWff(event.TLCat, f1, f2, p1, p2)

			if str(event.SName) not in nt2_npass :
				nt2_npass[str(event.SName)] = {}
				nt10_npass[str(event.SName)] = {}
				nt01_npass[str(event.SName)] = {}
				nt0_npass[str(event.SName)] = {}
				for ch_str, charge in self.charges.iteritems() :
					nt2_npass[str(event.SName)][ch_str] = {}
					nt10_npass[str(event.SName)][ch_str] = {}
					nt01_npass[str(event.SName)][ch_str] = {}
					nt0_npass[str(event.SName)][ch_str] = {}
					for chan in res[ch_str] :
						nt2_npass[str(event.SName)][ch_str][chan] = 0
						nt10_npass[str(event.SName)][ch_str][chan] = 0
						nt01_npass[str(event.SName)][ch_str][chan] = 0
						nt0_npass[str(event.SName)][ch_str][chan] = 0

			for ch_str, charge in self.charges.iteritems() :
				if charge != 0 and event.Charge != charge : continue

				# int
				res[ch_str]['al'].npp += npp;
				res[ch_str]['al'].npf += npf;
				res[ch_str]['al'].nfp += nfp;
				res[ch_str]['al'].nff += nff;
				if event.TLCat is 0 : res[ch_str]['al'].nt2  += scale;  nt2_npass[str(event.SName)][ch_str]['al'] += 1;
				if event.TLCat is 1 : res[ch_str]['al'].nt10 += scale; nt10_npass[str(event.SName)][ch_str]['al'] += 1;
				if event.TLCat is 2 : res[ch_str]['al'].nt01 += scale; nt01_npass[str(event.SName)][ch_str]['al'] += 1;
				if event.TLCat is 3 : res[ch_str]['al'].nt0  += scale;  nt0_npass[str(event.SName)][ch_str]['al'] += 1;

				# MM
				if event.Flavor is 0 :
					res[ch_str]['mm'].npp += npp;
					res[ch_str]['mm'].npf += npf;
					res[ch_str]['mm'].nfp += nfp;
					res[ch_str]['mm'].nff += nff;
					if event.TLCat is 0 : res[ch_str]['mm'].nt2  += scale;  nt2_npass[str(event.SName)][ch_str]['mm'] += 1;
					if event.TLCat is 1 : res[ch_str]['mm'].nt10 += scale; nt10_npass[str(event.SName)][ch_str]['mm'] += 1;
					if event.TLCat is 2 : res[ch_str]['mm'].nt01 += scale; nt01_npass[str(event.SName)][ch_str]['mm'] += 1;
					if event.TLCat is 3 : res[ch_str]['mm'].nt0  += scale;  nt0_npass[str(event.SName)][ch_str]['mm'] += 1;

				# EM
				if event.Flavor is 1 :
					res[ch_str]['em'].npp += npp
					res[ch_str]['em'].npf += npf
					res[ch_str]['em'].nfp += nfp
					res[ch_str]['em'].nff += nff
					if event.TLCat is 0 : res[ch_str]['em'].nt2  += scale;  nt2_npass[str(event.SName)][ch_str]['em'] += 1;
					if event.TLCat is 1 : res[ch_str]['em'].nt10 += scale; nt10_npass[str(event.SName)][ch_str]['em'] += 1;
					if event.TLCat is 2 : res[ch_str]['em'].nt01 += scale; nt01_npass[str(event.SName)][ch_str]['em'] += 1;
					if event.TLCat is 3 : res[ch_str]['em'].nt0  += scale;  nt0_npass[str(event.SName)][ch_str]['em'] += 1;

				# EE
				if event.Flavor is 2 :
					res[ch_str]['ee'].npp += npp
					res[ch_str]['ee'].npf += npf
					res[ch_str]['ee'].nfp += nfp
					res[ch_str]['ee'].nff += nff
					if event.TLCat is 0 : res[ch_str]['ee'].nt2  += scale;  nt2_npass[str(event.SName)][ch_str]['ee'] += 1;
					if event.TLCat is 1 : res[ch_str]['ee'].nt10 += scale; nt10_npass[str(event.SName)][ch_str]['ee'] += 1;
					if event.TLCat is 2 : res[ch_str]['ee'].nt01 += scale; nt01_npass[str(event.SName)][ch_str]['ee'] += 1;
					if event.TLCat is 3 : res[ch_str]['ee'].nt0  += scale;  nt0_npass[str(event.SName)][ch_str]['ee'] += 1;

		#####################
		# Fakes predictions #
		#####################

		print '[status] calculating fake predictions..'

		FR.setNToyMCs(100)
		FR.setAddESyst(self.FakeESyst)

		# ratios with ewk subtraction
		FR.setMFRatio(self.fpr.MufRatio_MC, self.fpr.MufRatioE_MC) # set error to pure statistical of ratio
		FR.setEFRatio(self.fpr.ElfRatio_MC, self.fpr.ElfRatioE_MC)
		FR.setMPRatio(self.fpr.MupRatio_MC, self.fpr.MupRatioE_MC)
		FR.setEPRatio(self.fpr.ElpRatio_MC, self.fpr.ElpRatioE_MC)

		for ch_str in self.charges :

			# store syst errors
			res[ch_str]['mm'].npp_systerr = self.FakeESyst*res[ch_str]['mm'].npp;
			res[ch_str]['em'].npp_systerr = self.FakeESyst*res[ch_str]['em'].npp;
			res[ch_str]['ee'].npp_systerr = self.FakeESyst*res[ch_str]['ee'].npp;

			res[ch_str]['mm'].npf_systerr = self.FakeESyst*res[ch_str]['mm'].npf;
			res[ch_str]['em'].npf_systerr = self.FakeESyst*res[ch_str]['em'].npf;
			res[ch_str]['ee'].npf_systerr = self.FakeESyst*res[ch_str]['ee'].npf;

			res[ch_str]['mm'].nfp_systerr = self.FakeESyst*res[ch_str]['mm'].nfp;
			res[ch_str]['em'].nfp_systerr = self.FakeESyst*res[ch_str]['em'].nfp;
			res[ch_str]['ee'].nfp_systerr = self.FakeESyst*res[ch_str]['ee'].nfp;

			res[ch_str]['mm'].nff_systerr = self.FakeESyst*res[ch_str]['mm'].nff;
			res[ch_str]['em'].nff_systerr = self.FakeESyst*res[ch_str]['em'].nff;
			res[ch_str]['ee'].nff_systerr = self.FakeESyst*res[ch_str]['ee'].nff;

			for s in nt2_npass :

				scale = self.lumi / self.samples[s].getLumi()
#				scale = 1.

				print '%5.2f %5.2f %5.2f %5.2f' % (nt2_npass[s][ch_str]['mm'], nt10_npass[s][ch_str]['mm'], nt01_npass[s][ch_str]['mm'], nt0_npass[s][ch_str]['mm'])
				print '%5.2f %5.2f %5.2f %5.2f' % (nt2_npass[s][ch_str]['em'], nt10_npass[s][ch_str]['em'], nt01_npass[s][ch_str]['em'], nt0_npass[s][ch_str]['em'])
				print '%5.2f %5.2f %5.2f %5.2f' % (nt2_npass[s][ch_str]['ee'], nt10_npass[s][ch_str]['ee'], nt01_npass[s][ch_str]['ee'], nt0_npass[s][ch_str]['ee'])

				print '%5.2f %5.2f %5.2f %5.2f' % (res[ch_str]['mm'].nt2, res[ch_str]['mm'].nt10, res[ch_str]['mm'].nt01, res[ch_str]['mm'].nt0)
				print '%5.2f %5.2f %5.2f %5.2f' % (res[ch_str]['em'].nt2, res[ch_str]['em'].nt10, res[ch_str]['em'].nt01, res[ch_str]['em'].nt0)
				print '%5.2f %5.2f %5.2f %5.2f' % (res[ch_str]['ee'].nt2, res[ch_str]['ee'].nt10, res[ch_str]['ee'].nt01, res[ch_str]['ee'].nt0)

				FR.setMMNtl(nt2_npass[s][ch_str]['mm'], nt10_npass[s][ch_str]['mm'], nt01_npass[s][ch_str]['mm'], nt0_npass[s][ch_str]['mm'])
				FR.setEMNtl(nt2_npass[s][ch_str]['em'], nt10_npass[s][ch_str]['em'], nt01_npass[s][ch_str]['em'], nt0_npass[s][ch_str]['em'])
				FR.setEENtl(nt2_npass[s][ch_str]['ee'], nt10_npass[s][ch_str]['ee'], nt01_npass[s][ch_str]['ee'], nt0_npass[s][ch_str]['ee'])

				# store stat errors
				res[ch_str]['mm'].npp_staterr = math.sqrt(res[ch_str]['mm'].npp_staterr**2 + (scale*FR.getMMNppEStat())**2);
				res[ch_str]['em'].npp_staterr = math.sqrt(res[ch_str]['em'].npp_staterr**2 + (scale*FR.getEMNppEStat())**2);
				res[ch_str]['ee'].npp_staterr = math.sqrt(res[ch_str]['ee'].npp_staterr**2 + (scale*FR.getEENppEStat())**2);

				res[ch_str]['mm'].npf_staterr = math.sqrt(res[ch_str]['mm'].npf_staterr**2 + (scale*FR.getMMNpfEStat())**2);
				res[ch_str]['em'].npf_staterr = math.sqrt(res[ch_str]['em'].npf_staterr**2 + (scale*FR.getEMNpfEStat())**2);
				res[ch_str]['ee'].npf_staterr = math.sqrt(res[ch_str]['ee'].npf_staterr**2 + (scale*FR.getEENpfEStat())**2);

				res[ch_str]['mm'].nfp_staterr = math.sqrt(res[ch_str]['mm'].nfp_staterr**2 + (scale*FR.getMMNfpEStat())**2);
				res[ch_str]['em'].nfp_staterr = math.sqrt(res[ch_str]['em'].nfp_staterr**2 + (scale*FR.getEMNfpEStat())**2);
				res[ch_str]['ee'].nfp_staterr = math.sqrt(res[ch_str]['ee'].nfp_staterr**2 + (scale*FR.getEENfpEStat())**2);

				res[ch_str]['mm'].nff_staterr = math.sqrt(res[ch_str]['mm'].nff_staterr**2 + (scale*FR.getMMNffEStat())**2);
				res[ch_str]['em'].nff_staterr = math.sqrt(res[ch_str]['em'].nff_staterr**2 + (scale*FR.getEMNffEStat())**2);
				res[ch_str]['ee'].nff_staterr = math.sqrt(res[ch_str]['ee'].nff_staterr**2 + (scale*FR.getEENffEStat())**2);

				res[ch_str]['al'].fake_staterr = math.sqrt(res[ch_str]['al'].fake_staterr**2 + (scale*FR.getTotEStat()  )**2)
				res[ch_str]['mm'].fake_staterr = math.sqrt(res[ch_str]['mm'].fake_staterr**2 + (scale*FR.getMMTotEStat())**2)
				res[ch_str]['em'].fake_staterr = math.sqrt(res[ch_str]['em'].fake_staterr**2 + (scale*FR.getEMTotEStat())**2)
				res[ch_str]['ee'].fake_staterr = math.sqrt(res[ch_str]['ee'].fake_staterr**2 + (scale*FR.getEETotEStat())**2)

			# store fake predictions
			res[ch_str]['al'].fake_err = math.sqrt(res[ch_str]['al'].fake_staterr**2 + self.FakeESyst2*res[ch_str]['al'].fake**2)
			res[ch_str]['mm'].fake_err = math.sqrt(res[ch_str]['mm'].fake_staterr**2 + self.FakeESyst2*res[ch_str]['mm'].fake**2)
			res[ch_str]['em'].fake_err = math.sqrt(res[ch_str]['em'].fake_staterr**2 + self.FakeESyst2*res[ch_str]['em'].fake**2)
			res[ch_str]['ee'].fake_err = math.sqrt(res[ch_str]['ee'].fake_staterr**2 + self.FakeESyst2*res[ch_str]['ee'].fake**2)

			# store tight-tight yields and errors
			for chan in res[ch_str] :
				staterr2 = 0.
				for s in nt2_npass :
					scale = self.lumi / self.samples[s].getLumi()
					#scale = 1.
					staterr = scale * self.samples[s].getError(nt2_npass[s][ch_str][chan])
					staterr2 += staterr**2
				res[ch_str][chan].nt2_staterr = math.sqrt(staterr2)

		return res


if __name__ == '__main__' :
	args = sys.argv
	selfile  = ''
	IntPred  = False
	DiffPred = False
	IntMC    = False
	DiffMC   = False
	RatioPlots        = False
	RatioControlPlots = False
	FakeClosure = False
	SLPlots = False

	if ('--help' in args) or ('-h' in args) or ('-d' not in args) :
		print 'usage: plotter.py -b <OPTIONS>'
		print ''
		print '\t-d <DIRECTORY OF SSDLYields.root>'
		print '\t-c <DATACARDFILE>'
		print '\t-s <SELECTIONFILE>'
		print ''
		print '\t--IntPred'
		print '\t\t- predicts integrated yields'
		print '\t\t- systematics'
		print '\t\t- creates datacards'
		print '\t\t- creates some tables'
		print ''
		print '\t--DiffPred'
		print '\t\t- predictions in various selections'
		print '\t\t- creates differential predictions for various varibles'
		print ''
		print '\t--IntMC'
		print '\t\t- uses MC predictions only'
		print '\t\t- creates yields table for various selections'
		print ''
		print '\t--DiffMC'
		print '\t\t- uses MC predictions only'
		print '\t\t- creates differential plots for various selections'
		print ''
		print '\t--RatioPlots'
		print '\t\tcreates prompt and fake ratio plots'
		print ''
		print '\t--RaioControlPlots'
		print '\t\tcreates control plots of various variables in QCD enriched region'
		print ''
		print '\t--FakeClosure'
		print '\t\truns closure test for fake lepton background'
		sys.exit(1)

	if ('-d' in args) and (args[args.index('-d')+1] != '') :
		path = str(args[args.index('-d')+1])
		print path

	if ('-c' in args) and (args[args.index('-c')+1] != '') :
		cardfile = str(args[args.index('-c')+1])
		print cardfile
	else :
		cardfile = '%s/DataCard_SSDL.dat' % path

	if ('-s' in args) and (args[args.index('-s')+1] != '') :
		selfile = str(args[args.index('-s')+1])
		print selfile

	if ('--IntPred' in args) :
		IntPred = True

	if ('--DiffPred' in args) :
		DiffPred = True

	if ('--IntMC' in args) :
		IntMC = True

	if ('--DiffMC' in args) :
		DiffMC = True

	if ('--RatioPlots' in args) :
		RatioPlots = True

	if ('--RatioControlPlots' in args) :
		if (len(args) > args.index('--RatioControlPlots')+1) and (not args[args.index('--RatioControlPlots')+1].startswith('-')) :
			RatioControlPlots = args[args.index('--RatioControlPlots')+1]
		else :
			print '[ERROR] Please specify which ratio control plots you like!'
			sys.exit(1)

	if ('--FakeClosure' in args) :
		FakeClosure = True

	if ('--SLPlots' in args) :
		SLPlots = True

	pl = plotter(path, cardfile, selfile)
	pl.do_analysis(IntPred = IntPred, DiffPred = DiffPred, IntMC = IntMC, DiffMC = DiffMC, RatioPlots = RatioPlots, RatioControlPlots = RatioControlPlots, FakeClosure = FakeClosure, SLPlots = SLPlots)
