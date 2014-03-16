#! /usr/bin/python

class selection :
	'''define and store selections'''

	def __init__(self, name, minHT = 0., maxHT = 8000., minMET = 0., maxMET = 8000, minNjets = 0, maxNjets = 99, minNbjetsL = 0, maxNbjetsL = 99, minNbjetsM = 0, maxNbjetsM = 99, minPt1 = 20., minPt2 = 20., applyZVeto = True, charge = 0, ttw = True, systflag = 0, flavor = -1) :
		self.name       = name
		self.minHT      = minHT
		self.maxHT      = maxHT
		self.minMET     = minMET
		self.maxMET     = maxMET
		self.minNjets   = minNjets
		self.maxNjets   = maxNjets
		self.minNbjetsL = minNbjetsL
		self.maxNbjetsL = maxNbjetsL
		self.minNbjetsM = minNbjetsM
		self.maxNbjetsM = maxNbjetsM
		self.minPt1     = minPt1
		self.minPt2     = minPt2
		self.applyZVeto = applyZVeto
		self.charge     = charge     # 0: all, +1: ++, -1: --
		self.ttw        = ttw
		self.systflag   = systflag
		self.flavor     = flavor     # -2: all, -1: same-sign, 0: mu-mu, 1: el-mu, 2: el-el, 3: mu-mu OS, 4: el-mu OS, 5: el-el OS

	def passes_selection(self, event, ttLeptons = True) :
		if event.HT     < self.minHT      : return False
		if event.HT     > self.maxHT      : return False
		if event.MET    < self.minMET     : return False
		if event.MET    > self.maxMET     : return False
		if event.NJ     < self.minNjets   : return False
		if event.NJ     > self.maxNjets   : return False
		if event.NbJ    < self.minNbjetsL : return False
		if event.NbJ    > self.maxNbjetsL : return False
		if event.NbJmed < self.minNbjetsM : return False
		if event.NbJmed > self.maxNbjetsM : return False
		if max(event.pT1, event.pT2) < self.minPt1  : return False
		if min(event.pT1, event.pT2) < self.minPt2  : return False
		if self.applyZVeto and event.PassZVeto is 0 : return False
		if self.charge is not 0 and event.Charge is not self.charge : return False
		if ttLeptons and event.TLCat > 0  : return False
		if self.flavor > -1  and event.Flavor is not self.flavor : return False
		if self.flavor is -1 and event.Flavor > 2                : return False
		if event.SystFlag is not self.systflag : return False
		return True

	def get_selectionString(self, ttLeptons = True) :
		'''return a selection string which can be used to draw directly from the tree'''
		selectionString = ''
		selectionString +=     'HT >= %f'     % (self.minHT     )
		selectionString += ' && HT <= %f'     % (self.maxHT     )
		selectionString += ' && MET >= %f'    % (self.minMET    )
		selectionString += ' && MET <= %f'    % (self.maxMET    )
		selectionString += ' && NJ >= %d'     % (self.minNjets  )
		selectionString += ' && NJ <= %d'     % (self.maxNjets  )
		selectionString += ' && NbJ >= %d'    % (self.minNbjetsL)
		selectionString += ' && NbJ <= %d'    % (self.maxNbjetsL)
		selectionString += ' && NbJmed >= %d' % (self.minNbjetsM)
		selectionString += ' && NbJmed <= %d' % (self.maxNbjetsM)
		selectionString += ' && TMath::Max(pT1,pT2) >= %f' % (self.minPt1)
		selectionString += ' && TMath::Min(pT1,pT2) >= %f' % (self.minPt2)
		selectionString += ' && SystFlag == %d'            % (self.systflag)
		if self.applyZVeto      : selectionString += ' && PassZVeto != 0'
		if self.charge is not 0 : selectionString += ' && Charge == %d' % (self.charge)
		if ttLeptons            : selectionString += ' && TLCat == 0'
		if self.flavor > -1     : selectionString += ' && Flavor == %d' % (self.flavor)
		if self.flavor is -1    : selectionString += ' && Flavor < 3'
		return selectionString
