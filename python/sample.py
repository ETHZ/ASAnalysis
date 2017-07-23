#! /usr/bin/python

import ROOT

class sample :
	'''stores informations and numbers about a sample'''


	def __init__(self, name, datamc = -1, channel = -1, xsec = -1, ngen = -1, inputfile = '') :
		self.name    = name
		self.datamc  = datamc
		self.channel = channel
		self.xsec    = xsec
		self.ngen    = ngen
		self.inputfile = inputfile

		self.qcd     = ['QCD', 'MuEnr15', 'MuEnr10', 'MuEnr20', 'MuEnr30']
		self.top     = ['SingleT', 'TTJets']
		self.ewk     = ['DYJets', 'GJets', 'WJets', 'WbbJets']
		self.rare    = ['HWW', 'HZZ', 'HTauTau', 'TTbarW', 'TTbarZ', 'TTbarG', 'TbZ', 'DPSWW', 'WWZ', 'WZZ', 'WZZ', 'WWG', 'ZZZ', 'WWW', 'W+W+', 'W-W-', 'TTbarWW']
		self.diboson = ['GVJets', 'WGstarMu', 'WGstarTau', 'WWTo2L2Nu', 'WZTo3LNu', 'ZZTo4L']

		# the order matters here:
		self.setType()
		self.setSampleType()


	def __str__(self) :
		return '%14s:\tdatamc: %2d, channel: %2d, xsec: %9.2e, ngen: %9d, lumi: %9.2e' % (self.name, self.datamc, self.channel, self.xsec, self.ngen, self.getLumi())


	def getLumi(self) :
		'''returns luminosity'''
		if self.datamc is 0 : return -1.
		if self.ngen > 0 and self.xsec > 0 : return float(self.ngen)/self.xsec
		else : return -1.


	def setType(self) :
		'''-1: undef, 0: data, 1: QCD, 2: top, 3: EWK, 4: Rare SM, 5: diboson'''
		if self.datamc == 0                              : self.type = 0
		elif any([i in self.name for i in self.qcd    ]) : self.type = 1
		elif any([i in self.name for i in self.top    ]) : self.type = 2
		elif any([i in self.name for i in self.ewk    ]) : self.type = 3
		elif any([i in self.name for i in self.rare   ]) : self.type = 4
		elif any([i in self.name for i in self.diboson]) : self.type = 5
		else :
			print '[error] type of %s is not defined!' % (self.name)
			self.type = -1


	def getType(self) :
		'''-1: undef, 0: data, 1: QCD, 2: top, 3: EWK, 4: Rare SM, 5: diboson'''
		if self.type == -1 :
			print '[error] type of %s is not defined!' % (self.name)
			return -1
		return self.type


	def setSampleType(self) :
		'''set sample type'''
		if   self.name.startswith('DoubleMu' ) : self.sampletype = 0
		elif self.name.startswith('DoubleEle') : self.sampletype = 1
		elif self.name.startswith('MuEG'     ) : self.sampletype = 2
		elif self.name.startswith('MuHad'    ) : self.sampletype = 3
		elif self.name.startswith('EleHad'   ) : self.sampletype = 4
		elif self.name.startswith('SingleMu' ) : self.sampletype = 5
		elif self.datamc > 0 and self.getType() != 4 and self.getType() != 5 : self.sampletype = 10  # SM MC without RARE and di-BOSON
		elif self.getType() == 4 or self.getType() == 5                      : self.sampletype = 15  # RARE MC + di-BOSON
		elif self.datamc == 2                                                : self.sampletype = 20  # signal
		else :
			print '[error] sample type of %s is not defined!' % (self.name)
			self.sampletype = -1


	def getSampleType(self) :
		'''return sample type'''
		if self.sampletype == -1 :
			print '[error] sample type of %s is not defined!' % (self.name)
			return -1
		return self.sampletype


	def getError(self, n) :
		'''If n passed of ngen generated, what is upper limit on number of events passing?'''
		if self.ngen <= 0 : return -1.
		eff = ROOT.TEfficiency()
		upper = eff.ClopperPearson(self.ngen, int(n), 0.68, True)
		delta = upper - float(n)/float(self.ngen)
		return delta * float(self.ngen)


	def getError2(self, n) :
		err = self.getError(n)
		return err*err


	@staticmethod
	def get_samples(channel, samples) :
		include_DYJets10To50 = False
		if channel == 'DYJets' or channel == 'zjets' :
			if include_DYJets10To50 : return ['DYJets', 'DYJets10To50']
			else                    : return ['DYJets']
		if channel == 'WJets'  or channel == 'wjets' : return ['WJets' ]
		if channel == 'TTW'    or channel == 'ttw' : return ['TTbarW']
		if channel == 'TTZ'    or channel == 'ttz' : return ['TTbarZ']
		if channel == 'WZ'     or channel == 'wz' : return ['WZTo3LNu']
		if channel == 'TTH'    or channel == 'tth' : return ['HWW', 'HZZ', 'HTauTau']
		if channel == 'TTJets' or channel == 'ttbar'  : return [s for s in samples if s.startswith('TTJets')]
		if channel == 'Multiboson' : return ['WWW', 'WWZ', 'WZZ', 'ZZZ', 'WWG']
		samplelist = []
		for name, sample in samples.iteritems() :
			if not include_DYJets10To50 and sample.name == 'DYJets10To50' : continue
			if   channel == 'DoubleMu' :
				if (sample.datamc == 0) and (sample.channel == 0) : samplelist.append(sample.name)
			elif channel == 'DoubleEle' :
				if (sample.datamc == 0) and (sample.channel == 1) : samplelist.append(sample.name)
			elif channel == 'MuEG' :
				if (sample.datamc == 0) and (sample.channel == 2) : samplelist.append(sample.name)
			elif channel == 'SingleMu' :
				if (sample.datamc == 0) and (sample.channel == 5) : samplelist.append(sample.name)
			elif channel == 'SingleDoubleMu' :
				if (sample.datamc == 0) and ((sample.channel == 0) or (sample.channel == 5)) : samplelist.append(sample.name)
			elif channel == 'MC' :
				if (sample.datamc > 0) : samplelist.append(sample.name)
			elif channel == 'QCD' or channel == 'qcd' :
				if (sample.getType() == 1) : samplelist.append(sample.name)
			elif channel == 'Top' or channel == 'top' :
				if (sample.getType() == 2) : samplelist.append(sample.name)
			elif channel == 'Rare' or channel == 'rare' :
				if (sample.getSampleType() == 15) and (sample.name != 'TTbarW') and (sample.name != 'TTbarZ') and (sample.name != 'WZTo3LNu') : samplelist.append(sample.name)
		return samplelist
