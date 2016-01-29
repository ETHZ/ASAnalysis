#! /usr/bin/python
import sys
import ROOT
import plotter
import helper
import ttvStyle
import sample
import numpy


class analyzer :
	'''runs over analysis tree'''


	def __init__(self, path, ssdl_path, cardfile, suffix = '') :
		if not path.endswith('/') :
			path += '/'
		self.path = path
		helper.mkdir(self.path)
		self.pl = plotter.plotter(ssdl_path, cardfile)
		self.samples = self.pl.readDatacard(cardfile = cardfile, verbose = 0, event_count = True)
		for s in self.samples :
			print self.samples[s]
		self.histos = {}
		if not suffix.startswith('_') and suffix != '' :
			suffix = '_' + suffix
		file_name = 'SLYields%s.root' % suffix
		print '[status] generating output file %s..' % file_name
		self.slFile = ROOT.TFile('%s%s' % (self.path, file_name), 'RECREATE')
		self.slTree = ROOT.TTree('SLEvents', 'SLEventsTree')
		self.book_slTree()


	def isLooseMuon(self, event, mu) :
		if mu >= event.NMus : return False
		if event.MuPt[mu] < 20. : return False
		if event.MuEMVetoEt[mu]  > 4.0 : return False
		if event.MuHadVetoEt[mu] > 6.0 : return False
		if event.MuPassesTightID[mu] != 1 : return False
		if event.MuPFIso[mu] > 1.00 : return False
		if abs(event.MuD0[mu]) > 0.005 : return False
		return True


	def isTightMuon(self, event,  mu) :
		if not self.isLooseMuon(event, mu) : return False
		if event.MuPFIso[mu] > 0.05 : return False
		return True


	def isLooseElectron(self, event, el) :
		if event.ElPt[el] < 20. : return False
		for mu in range(event.NMus) :
			if not self.isLooseMuon(event, mu) : continue
			if not self.isTightMuon(event, mu) : continue
			if helper.get_deltaR(event.MuEta[mu], event.ElEta[el], event.MuPhi[mu], event.ElPhi[el]) > 0.1 : continue
			return False
		if event.ElChIsCons[el] != 1 : return False
		if abs(event.ElSCEta[el]) > 1.4442 and abs(event.ElSCEta[el]) < 1.566 : return False
		if event.ElPFIso[el] > 0.6 : return False
		if abs(event.ElD0[el]) > 0.01 : return False
		return True


	def isTightElectron(self, event, el) :
		if not self.isLooseElectron(event, el) : return False
		if event.ElPFIso[el] > 0.05 : return False
		if event.ElIsGoodElId_MediumWP[el] != 1 : return False
		return True


	def get_looseMuonIndices(self, event) :
		looseMus = []
		for mu in range(event.NMus) :
			if self.isLooseMuon(event, mu) : looseMus.append(mu)
		return looseMus


	def get_looseElectronIndices(self, event) :
		looseEls = []
		for el in range(event.NEls) :
			if self.isLooseElectron(event, el) : looseEls.append(el)
		return looseEls


	def do_analysis(self, s, file_name = '') :
#		for s in self.samples :
#			self.analyze_tree(s)
		self.analyze_tree(s, file_name)
		self.slFile.Close()


	def analyze_tree(self, s, file_name = '') :
		self.book_histos(s)
		if self.samples[s].inputfile.endswith('.root') :
			inputfile = self.samples[s].inputfile
		elif file_name != '' :
			if not file_name.endswith('.root') :
				file_name += '.root'
			inputfile = self.samples[s].inputfile + file_name
		else :
			print '[ERROR] Check input file!'
			sys.exit(1)
		file = ROOT.TFile.Open(inputfile, 'READ')
		tree = file.Get('Analysis')

		self.slFile.cd()

		print '[status] looping over Analysis tree with %d events..' % tree.GetEntries()
		for index, event in enumerate(tree) :
			looseMus = self.get_looseMuonIndices(event)
			looseEls = self.get_looseElectronIndices(event)
			if len(looseMus) == 0 and len(looseEls) == 0 : continue

#			self.slTree_Weight  [0] = event.Weight
			self.slTree_SName.assign(s)
			self.slTree_Run             [0] = event.Run
			self.slTree_LS              [0] = event.LumiSec
			self.slTree_Event           [0] = event.Event
			self.slTree_MET             [0] = event.pfMET
			self.slTree_NVrtx           [0] = event.NVrtx
			self.slTree_MuPt            [0] = -9999.
			self.slTree_MuPFIso         [0] = -9999.
			self.slTree_MuEta           [0] = -9999.
			self.slTree_MuPhi           [0] = -9999.
			self.slTree_MuCharge        [0] = 0
			self.slTree_MuD0            [0] = -9999.
			self.slTree_MuDz            [0] = -9999.
			self.slTree_NMus            [0] = len(looseMus)
			self.slTree_IsSignalMuon    [0] = 0
			self.slTree_ElPt            [0] = -9999.
			self.slTree_ElPFIso         [0] = -9999.
			self.slTree_ElEta           [0] = -9999.
			self.slTree_ElPhi           [0] = -9999.
			self.slTree_ElCharge        [0] = 0
			self.slTree_ElD0            [0] = -9999.
			self.slTree_ElDz            [0] = -9999.
			self.slTree_NEls            [0] = len(looseEls)
			self.slTree_IsSignalElectron[0] = 0

			if len(looseMus) > 0 :
				self.slTree_MuPt            [0] = event.MuPt            [looseMus[0]]
				self.slTree_MuPFIso         [0] = event.MuPFIso         [looseMus[0]]
				self.slTree_MuD0            [0] = event.MuD0            [looseMus[0]]
				self.slTree_IsSignalMuon    [0] = event.IsSignalMuon    [looseMus[0]]
				self.slTree_MuEta           [0] = event.MuEta           [looseMus[0]]
				self.slTree_MuPhi           [0] = event.MuPhi           [looseMus[0]]
				self.slTree_MuCharge        [0] = event.MuCharge        [looseMus[0]]
				self.slTree_MuDz            [0] = event.MuDz            [looseMus[0]]

			if len(looseEls) > 0 :
				self.slTree_ElPt            [0] = event.ElPt            [looseEls[0]]
				self.slTree_ElPFIso         [0] = event.ElPFIso         [looseEls[0]]
				self.slTree_ElD0            [0] = event.ElD0            [looseEls[0]]
				self.slTree_IsSignalElectron[0] = event.IsSignalElectron[looseEls[0]]
				self.slTree_ElEta           [0] = event.ElEta           [looseEls[0]]
				self.slTree_ElPhi           [0] = event.ElPhi           [looseEls[0]]
				self.slTree_ElCharge        [0] = event.ElCh            [looseEls[0]]
				self.slTree_ElDz            [0] = event.ElDz            [looseEls[0]]

			if len(looseMus) > 0 or len(looseEls) > 0 : self.slTree.Fill()

#			# muon isolation
#			if len(looseMus) == 1 :
#				self.histos[s]['MuIso'].Fill(event.MuPFIso[looseMus[0]])
#
#			# electron isolation
#			if len(looseEls) == 1 :
#				self.histos[s]['ElIso'].Fill(event.ElPFIso[looseEls[0]])

		self.slTree.Write()
#		self.save_histos(s)



	def book_histos(self, s) :
		self.histos[s] = {}
		self.histos[s]['MuIso'] = ROOT.TH1D('%s_MuIso' % s, '%s MuIso' % s, 20, 0., 1.); self.histos[s]['MuIso'].Sumw2()
		self.histos[s]['ElIso'] = ROOT.TH1D('%s_ElIso' % s, '%s ElIso' % s, 20, 0., 1.); self.histos[s]['ElIso'].Sumw2()


	def add_histos(self, process) :
		samples = sample.sample.get_samples(process, self.samples)
		self.book_histos(process)

		for s in samples :
			for histo_name in self.histos[s] :
				self.histos[process][histo_name].Add(self.histos[s][histo_name], self.pl.lumi / self.samples[s].getLumi())


	def save_histos(self, s, suffix = '') :
		histofile = ROOT.TFile.Open('%s%s_histos%s.root' % ('./', s, suffix), 'RECREATE')
		histofile.cd()
		histofile.mkdir('%s' % s)
		histofile.cd   ('%s' % s)
		for name, histo in self.histos[s].iteritems() :
			histo.Write()
		histofile.Close()


	def book_slTree(self) :
#		self.slTree_Weight    = numpy.zeros(1, dtype=float); self.slTree.Branch('Weight',      self.slTree_Weight     , 'Weight/D'   );
		self.slTree_SName     = ROOT.std.string(          ); self.slTree.Branch('SName'      , self.slTree_SName);
		self.slTree_Run       = numpy.zeros(1, dtype=int  ); self.slTree.Branch('Run'        , self.slTree_Run        , 'Run/I'      );
		self.slTree_LS        = numpy.zeros(1, dtype=int  ); self.slTree.Branch('LS'         , self.slTree_LS         , 'LS/I'       );
		self.slTree_Event     = numpy.zeros(1, dtype=int  ); self.slTree.Branch('Event'      , self.slTree_Event      , 'Event/I'    );
		self.slTree_MET       = numpy.zeros(1, dtype=float); self.slTree.Branch('MET'        , self.slTree_MET        , 'MET/D'      );
		self.slTree_MuPt      = numpy.zeros(1, dtype=float); self.slTree.Branch('MuPt'       , self.slTree_MuPt       , 'MuPt/D'     );
		self.slTree_ElPt      = numpy.zeros(1, dtype=float); self.slTree.Branch('ElPt'       , self.slTree_ElPt       , 'ElPt/D'     );
		self.slTree_MuPFIso   = numpy.zeros(1, dtype=float); self.slTree.Branch('MuPFIso'    , self.slTree_MuPFIso    , 'MuPFIso/D'  );
		self.slTree_ElPFIso   = numpy.zeros(1, dtype=float); self.slTree.Branch('ElPFIso'    , self.slTree_ElPFIso    , 'ElPFIso/D'  );
		self.slTree_MuD0      = numpy.zeros(1, dtype=float); self.slTree.Branch('MuD0'       , self.slTree_MuD0       , 'MuD0/D'     );
		self.slTree_ElD0      = numpy.zeros(1, dtype=float); self.slTree.Branch('ElD0'       , self.slTree_ElD0       , 'ElD0/D'     );
		self.slTree_NVrtx     = numpy.zeros(1, dtype=int  ); self.slTree.Branch('NVrtx'      , self.slTree_NVrtx      , 'NVrtx/I'    );
		self.slTree_NMus      = numpy.zeros(1, dtype=int  ); self.slTree.Branch('NMus'       , self.slTree_NMus       , 'NMus/I'     );
		self.slTree_NEls      = numpy.zeros(1, dtype=int  ); self.slTree.Branch('NEls'       , self.slTree_NEls       , 'NEls/I'     );
		self.slTree_IsSignalMuon     = numpy.zeros(1, dtype=int  ); self.slTree.Branch('IsSignalMuon'    , self.slTree_IsSignalMuon    , 'IsSignalMuon/I'    );
		self.slTree_IsSignalElectron = numpy.zeros(1, dtype=int  ); self.slTree.Branch('IsSignalElectron', self.slTree_IsSignalElectron, 'IsSignalElectron/I');

		self.slTree_MuEta     = numpy.zeros(1, dtype=float); self.slTree.Branch('MuEta'      , self.slTree_MuEta      , 'MuEta/D'    );
		self.slTree_ElEta     = numpy.zeros(1, dtype=float); self.slTree.Branch('ElEta'      , self.slTree_ElEta      , 'ElEta/D'    );
		self.slTree_MuPhi     = numpy.zeros(1, dtype=float); self.slTree.Branch('MuPhi'      , self.slTree_MuPhi      , 'MuPhi/D'    );
		self.slTree_ElPhi     = numpy.zeros(1, dtype=float); self.slTree.Branch('ElPhi'      , self.slTree_ElPhi      , 'ElPhi/D'    );
		self.slTree_MuCharge  = numpy.zeros(1, dtype=int  ); self.slTree.Branch('MuCharge'   , self.slTree_MuCharge   , 'MuCharge/I' );
		self.slTree_ElCharge  = numpy.zeros(1, dtype=int  ); self.slTree.Branch('ElCharge'   , self.slTree_ElCharge   , 'ElCharge/I' );
		self.slTree_MuDz      = numpy.zeros(1, dtype=float); self.slTree.Branch('MuDz'       , self.slTree_MuDz       , 'MuDz/D'     );
		self.slTree_ElDz      = numpy.zeros(1, dtype=float); self.slTree.Branch('ElDz'       , self.slTree_ElDz       , 'ElDz/D'     );


if __name__ == '__main__' :
	args = sys.argv

	if '--help' in args or '-d' not in args :
		print 'usage: analyzer.py -b <OPTIONS>'
		print ''
		print '\t-d <DIRECTORY OF SSDLYields.root>'
		print '\t-c <DATACARDFILE>'
		print '\t-o <OUTPUT DIRECTORY>'
		print '\t--suffix <SUFFIX FOR SLYields.root>'
		sys.exit(1)

	if ('-d' in args) and (args[args.index('-d')+1] != '') :
		path = str(args[args.index('-d')+1])
		if not path.endswith('/') :
			path += '/'
		print '[info] directory of SSDLYields.root: %s' % path

	if ('-c' in args) and (args[args.index('-c')+1] != '') :
		cardfile = str(args[args.index('-c')+1])
	else :
		cardfile = path + 'DataCard_SSDL.dat'
	print '[info] cardfile: %s' % cardfile

	if ('-o' in args) and (args[args.index('-o')+1] != '') :
		output_path = str(args[args.index('-o')+1])
		if not output_path.endswith('/') :
			output_path += '/'
		print '[info] output path: %s' % output_path


	if '--suffix' in args and args[args.index('--suffix')+1] != '' :
		suffix = str(args[args.index('--suffix')+1])
	else :
		suffix = ''

	if '--sample' in args and args[args.index('--sample')+1] != '' :
		sample = str(args[args.index('--sample')+1])
	else :
		sample = ''

	if '--filename' in args and args[args.index('--filename')+1] != '' :
		filename = str(args[args.index('--filename')+1])
		suffix = filename
	else :
		filename = ''

	dp = analyzer(path = output_path, ssdl_path = path, cardfile = cardfile, suffix = 'test')
	dp.do_analysis(sample, filename)
