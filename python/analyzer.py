#! /usr/bin/python
import ROOT
import plotter
import helper
import ttvStyle
import sample


class analyzer :
	'''runs over analysis tree'''


	def __init__(self, path, cardfile) :
		self.path = path
		self.pl = plotter.plotter(path, cardfile)
		self.samples = self.pl.readDatacard(cardfile = cardfile, verbose = 0, event_count = True)
		for s in self.samples :
			print self.samples[s], self.samples[s].ngen
		self.histos = {}


	def isLooseMuon(self, event, mu) :
		if mu >= event.NMus : return False
		if event.MuPt[mu] < 20. : return False
		if event.MuEMVetoEt[mu]  > 4.0 : return False
		if event.MuHadVetoEt[mu] > 6.0 : return False
		if event.MuPassesTightID[mu] != 1 : return False
		if event.MuPFIso[mu] > 1.00 : return False # using detector isolation for ttWZ as requested by f.p.
		if abs(event.MuD0[mu]) > 0.005 : return False # this is for testing only!!!
		return True


	def isTightMuon(self, event,  mu) :
		if not self.isLooseMuon(event, mu) : return False
		if event.MuPFIso[mu] > 0.05 : return False # using detector isolation for ttWZ as requested by f.p.
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


	def do_analysis(self) :
		for s in self.samples :
			self.analyze_tree(s)
		self.add_histos('qcd')
		pl = ttvStyle.ttvStyle(lumi = 0., TeX_switch = False)
		canvas = pl.get_canvas('iso')
		canvas.cd()
		self.histos['qcd']['MuIso'].Draw('PE')
		canvas.Print('iso.pdf')
		helper.save_object(self.histos, './histos.pkl')


	def analyze_tree(self, s) :
		self.book_histos(s)
		file = ROOT.TFile.Open(self.samples[s].inputfile, 'READ')
		tree = file.Get('Analysis')

		for index, event in enumerate(tree) :
			if index < 20 :
				if event.NEls != 0 :
					print event.NEls, event.ElPFIso[0], self.get_looseElectronIndices(event)

			# muon isolation
			looseMus = self.get_looseMuonIndices(event)
			if len(looseMus) == 1 :
				self.histos[s]['MuIso'].Fill(event.MuPFIso[looseMus[0]])

			# electron isolation
			looseEls = self.get_looseElectronIndices(event)
			if len(looseEls) == 1 :
				self.histos[s]['ElIso'].Fill(event.ElPFIso[looseEls[0]])

		self.save_histos(s)



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
		self.slTree_Weight    = np.zeros(1, dtype=float); self.slTree.Branch('Weight',      self.slTree_Weight     , 'Weight/D'   );
		self.slTree_SName     = ROOT.std.string(       ); self.slTree.Branch('SName',       self.slTree_SName);
		self.slTree_SType     = np.zeros(1, dtype=int  ); self.slTree.Branch('SType',       self.slTree_SType      , 'SType/I'    );
		self.slTree_Run       = np.zeros(1, dtype=int  ); self.slTree.Branch('Run',         self.slTree_Run        , 'Run/I'      );
		self.slTree_LS        = np.zeros(1, dtype=int  ); self.slTree.Branch('LS',          self.slTree_LS         , 'LS/I'       );
		self.slTree_Event     = np.zeros(1, dtype=int  ); self.slTree.Branch('Event',       self.slTree_Event      , 'Event/I'    );
		self.slTree_HT        = np.zeros(1, dtype=float); self.slTree.Branch('HT',          self.slTree_HT         , 'HT/D'       );
		self.slTree_MET       = np.zeros(1, dtype=float); self.slTree.Branch('MET',         self.slTree_MET        , 'MET/D'      );
		self.slTree_MuPt      = np.zeros(1, dtype=float); self.slTree.Branch('MuPt',        self.slTree_MuPt       , 'MuPt/D'     );
		self.slTree_ElPt      = np.zeros(1, dtype=float); self.slTree.Branch('ElPt',        self.slTree_ElPt       , 'ElPt/D'     );
		self.slTree_MuPFIso   = np.zeros(1, dtype=float); self.slTree.Branch('MuPFIso',     self.slTree_MuPFIso    , 'MuPFIso/D'  );
		self.slTree_ElPFIso   = np.zeros(1, dtype=float); self.slTree.Branch('ElPFIso',     self.slTree_ElPFIso    , 'ElPFIso/D'  );
#		self.slTree_D01       = np.zeros(1, dtype=float); self.slTree.Branch('D01',         self.slTree_D01        , 'D01/D'      );
#		self.slTree_D02       = np.zeros(1, dtype=float); self.slTree.Branch('D02',         self.slTree_D02        , 'D02/D'      );
#		self.slTree_Rho       = np.zeros(1, dtype=float); self.slTree.Branch('Rho',         self.slTree_Rho        , 'Rho/D'      );
		self.slTree_NVrtx     = np.zeros(1, dtype=int  ); self.slTree.Branch('NVrtx',       self.slTree_NVrtx      , 'NVrtx/I'    );
		self.slTree_NMus      = np.zeros(1, dtype=int  ); self.slTree.Branch('NMus',        self.slTree_NMus       , 'NMus/I'     );
		self.slTree_NEls      = np.zeros(1, dtype=int  ); self.slTree.Branch('NEls',        self.slTree_NEls       , 'NEls/I'     );


if __name__ == '__main__' :
	dp = analyzer('../plots/Oct11-DYJets10To50-SystStudies-looseOS-Apr16', '../DataCard_SSDL.dat')
	dp.do_analysis()
