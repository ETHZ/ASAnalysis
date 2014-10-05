#! /usr/bin/python

class selection :
	'''define and store selections'''


	def __init__(self, name, minHT = 0., maxHT = 8000., minMET = 0., maxMET = 8000, minNjets = 0, maxNjets = 99, minNbjetsL = 0, maxNbjetsL = 99, minNbjetsM = 0, maxNbjetsM = 99, minPt1 = 20., minPt2 = 20., minMTLep1 = 0., minMTLep2 = 0., maxMTLep1 = 8000., maxMTLep2 = 8000., applyZVeto = True, charge = 0, ttw = True, systflag = 0, flavor = -1, mll = 8., sname = '') :
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
		self.minMTLep1  = minMTLep1
		self.minMTLep2  = minMTLep2
		self.maxMTLep1  = maxMTLep1
		self.maxMTLep2  = maxMTLep2
		self.ZVeto      = applyZVeto
		self.charge     = charge     # 0: all, +1: ++, -1: --
		self.ttw        = ttw
		self.systflag   = systflag
		self.flavor     = flavor     # -2: all, -1: same-sign, 0: mu-mu, 1: el-mu, 2: el-el, 3: mu-mu OS, 4: el-mu OS, 5: el-el OS
		self.mll        = mll
		self.sname      = sname


	def __str__(self) :
		printout  = '\n%s applies the following cuts:\n' % (self.name)
		printout += '\n   %5.1f <= HT                <= %6.1f' % (self.minHT     , self.maxHT     )
		printout += '\n   %5.1f <= MET               <= %6.1f' % (self.minMET    , self.maxMET    )
		printout += '\n   %3d   <= N jets            <= %4d'   % (self.minNjets  , self.maxNjets  )
		printout += '\n   %3d   <= N b-tags loose    <= %4d'   % (self.minNbjetsL, self.maxNbjetsL)
		printout += '\n   %3d   <= N b-tags medium   <= %4d'   % (self.minNbjetsM, self.maxNbjetsM)
		printout += '\n   %5.1f <= lead. lept. pT'             % (self.minPt1    )
		printout += '\n   %5.1f <= sublead. lept. pT'          % (self.minPt2    )
		if self.ZVeto       : printout += '\n   Z veto'
		if self.charge > 0  : printout += '\n   only l+l+ events'
		if self.charge < 0  : printout += '\n   only l-l- events'
		if self.flavor > -1 : printout += '\n   flavor = %d' % (self.flavor)
		printout                       += '\n   systflag = %d' % (self.systflag)
		printout += '\n'
		return printout


	def set_cutsFromDict(self, cuts) :
		if 'name'       in cuts.keys() : self.name       = cuts['name'      ]
		if 'minHT'      in cuts.keys() : self.minHT      = cuts['minHT'     ]
		if 'maxHT'      in cuts.keys() : self.maxHT      = cuts['maxHT'     ]
		if 'minMET'     in cuts.keys() : self.minMET     = cuts['minMET'    ]
		if 'maxMET'     in cuts.keys() : self.maxMET     = cuts['maxMET'    ]
		if 'minNjets'   in cuts.keys() : self.minNjets   = cuts['minNjets'  ]
		if 'maxNjets'   in cuts.keys() : self.maxNjets   = cuts['maxNjets'  ]
		if 'minNbjetsL' in cuts.keys() : self.minNbjetsL = cuts['minNbjetsL']
		if 'maxNbjetsL' in cuts.keys() : self.maxNbjetsL = cuts['maxNbjetsL']
		if 'minNbjetsM' in cuts.keys() : self.minNbjetsM = cuts['minNbjetsM']
		if 'maxNbjetsM' in cuts.keys() : self.maxNbjetsM = cuts['maxNbjetsM']
		if 'minPt1'     in cuts.keys() : self.minPt1     = cuts['minPt1'    ]
		if 'minPt2'     in cuts.keys() : self.minPt2     = cuts['minPt2'    ]
		if 'minMTLep1'  in cuts.keys() : self.minMTLep1  = cuts['minMTLep1' ]
		if 'minMTLep2'  in cuts.keys() : self.minMTLep2  = cuts['minMTLep2' ]
		if 'maxMTLep1'  in cuts.keys() : self.maxMTLep1  = cuts['maxMTLep1' ]
		if 'maxMTLep2'  in cuts.keys() : self.maxMTLep2  = cuts['maxMTLep2' ]
		if 'applyZVeto' in cuts.keys() : self.ZVeto      = cuts['applyZVeto']
		if 'charge'     in cuts.keys() : self.charge     = cuts['charge'    ]
		if 'ttw'        in cuts.keys() : self.ttw        = cuts['ttw'       ]
		if 'systflag'   in cuts.keys() : self.systflag   = cuts['systflag'  ]
		if 'flavor'     in cuts.keys() : self.flavor     = cuts['flavor'    ]
		if 'mll'        in cuts.keys() : self.mll        = cuts['mll'       ]
		if 'sname'      in cuts.keys() : self.sname      = cuts['sname'     ]


	def passes_selection(self, event, ttLeptons = True, noChargeSel = False, OSwoZVeto = False) :
		if event.SystFlag != self.systflag                                              : return False
		if not (OSwoZVeto and event.Flavor > 2) and self.ZVeto and event.PassZVeto == 0 : return False
		if self.charge != 0 and event.Charge != self.charge and not (noChargeSel)       : return False
		if ttLeptons and event.TLCat > 0                                                : return False
		if self.flavor > -1  and event.Flavor != self.flavor                            : return False
		if self.flavor == -1 and event.Flavor > 2 and not (OSwoZVeto)                   : return False
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
		if event.Mll    < self.mll        : return False
		if max(event.pT1, event.pT2) < self.minPt1 : return False
		if min(event.pT1, event.pT2) < self.minPt2 : return False
		if event.MTLep1 < self.minMTLep1 : return False
		if event.MTLep2 < self.minMTLep2 : return False
		if event.MTLep1 > self.maxMTLep1 : return False
		if event.MTLep2 > self.maxMTLep2 : return False
		if self.sname != '' and self.sname != str(event.SName) : return False

		return True


	def get_selectionString(self, OS_data = (-1, -1), ttLeptons = True, ResTree = False) :
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
		selectionString += ' && Mll >= %f'    % (self.mll       )
		selectionString += ' && TMath::Max(pT1,pT2) >= %f' % (self.minPt1)
		selectionString += ' && TMath::Min(pT1,pT2) >= %f' % (self.minPt2)
		selectionString += ' && MTLep1 >= %f' % (self.minMTLep1 )
		selectionString += ' && MTLep2 >= %f' % (self.minMTLep2 )
		selectionString += ' && MTLep1 <= %f' % (self.maxMTLep1 )
		selectionString += ' && MTLep2 <= %f' % (self.maxMTLep2 )
		selectionString += ' && SystFlag == %d'            % (self.systflag)
		if not ResTree :
			if OS_data[0] < 0 :
				if self.ZVeto        : selectionString += ' && PassZVeto != 0'
				if self.charge != 0  : selectionString += ' && Charge == %d' % (self.charge)
				if ttLeptons         : selectionString += ' && TLCat == 0'
				if self.flavor > -1  : selectionString += ' && Flavor == %d' % (self.flavor)
				if self.flavor == -1 : selectionString += ' && Flavor < 3'
			else :
				if OS_data[1] < 0 : return -1
				selectionString += ' && SType < 3'
				selectionString += ' && Flavor == %d' % (OS_data[0])
				if ttLeptons         : selectionString += ' && TLCat == 0'
				if OS_data[1] == 4 :
					selectionString += ' && (BECat == 1 || BECat == 2)'
				else :
					selectionString += ' && BECat == %d'  % (OS_data[1])
				if self.flavor > -1 : selectionString += ' && Flavor == %d' % (self.flavor+3)
		else :
			if self.flavor > -1                    : selectionString += ' && Flavor == %d' % (self.flavor)
			if OS_data[0] < 0 and self.charge != 0 : selectionString += ' && Charge == %d' % (self.charge)
		if self.sname != '' : selectionString += ' && SName == \"%s\"' % (self.sname)

		return selectionString
