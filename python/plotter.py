#! /usr/bin/python
#import os, sys, commands, subprocess, math
import ROOT
import sys
import math
import selection as sel
import sample
import helper
import ratios
import prediction
import result
import time


class plotter :
	'''the plotter reads sigtree and produces plots'''

	def __init__(self, path) :
		print '[status] initialize plotter..'
		self.path = path
		self.ssdlfile = ROOT.TFile.Open(path + '/SSDLYields.root', 'READ')
		if self.ssdlfile == None :
			sys.exit(1)
		self.sigtree = self.ssdlfile.Get('SigEvents')
		print '[status] loaded SigEventsTree with %d events' % (self.sigtree.GetEntries())

		self.skimfile = ROOT.TFile.Open(path + '/SSDLYields_skim.root', 'READ')
		if self.skimfile == None :
			print '[status] creating skimmed tree file..'
			self.skimfile = ROOT.TFile.Open(path + '/SSDLYields_skim.root', 'RECREATE')
			self.skimfile.cd()
			self.skimtree = self.sigtree.CopyTree('Flavor < 3 && (SType < 3 || TLCat == 0) && (SType > 2 || SystFlag == 0)') # Only same-sign, if MC only tight-tight. Only MC for syst studies.
			self.skimtree.AutoSave()
			self.skimfile.Write()
		else :
			self.skimtree = self.skimfile.Get('SigEvents')

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
		presel = sel.selection(name = 'presel', minNjets = 3, minNbjetsL = 1, minNbjetsM = 1)
		#print presel.get_selectionString()
		self.selections[presel.name] = presel

		self.chmid_sf = 1.62

		EWK_SF = {}
		EWK_SF['el']   = self.get_EWK_SF('el')
		EWK_SF['mu17'] = self.get_EWK_SF('mu17')
		EWK_SF['mu24'] = self.get_EWK_SF('mu24')

		applyEwkSubtr = True
#		self.fill_ratios(self.get_samples('SingleDoubleMu'), self.get_samples('DoubleEle'), 0, True, EWK_SF)
		self.fpr = ratios.ratios(path, self.samples)
		self.fpr.fill_ratios(self.get_samples('SingleDoubleMu'), self.get_samples('DoubleEle'), 0, True, EWK_SF)

		print 'MufRatio: %f +/- %f' % (self.fpr.MufRatio, self.fpr.MufRatioE)
		print 'ElfRatio: %f +/- %f' % (self.fpr.ElfRatio, self.fpr.ElfRatioE)
		print 'MupRatio: %f +/- %f' % (self.fpr.MupRatio, self.fpr.MupRatioE)
		print 'ElpRatio: %f +/- %f' % (self.fpr.ElpRatio, self.fpr.ElpRatioE)

		print "self.get_fRatio('Muon', 35., 0.5, 0):", self.fpr.get_fRatio('Muon', 35., 0.5, 0)
		print "self.get_fRatio('Muon', 55., 0.5, 0):", self.fpr.get_fRatio('Muon', 55., 0.5, 0)
		print "self.get_pRatio('Muon', 55., 0     ):", self.fpr.get_pRatio('Muon', 55., 0     )

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
#		res_syst = {}
#		res_syst['Normal'] = self.make_IntPredictions(presel)
#
#		return

#		sels = {}
#		systres = {}
#		for syst in self.systematics :
#			sels[syst] = sel.selection(name = 'final_' + syst, minNjets = 3, minNbjetsL = 1, minNbjetsM = 1, minPt1 = 40., minPt2 = 40., minHT = 155., systflag = self.systematics[syst])
#			systres[syst] = self.make_IntPredictions(sels[syst])
#
#		helper.save_obj(systres, 'systs')
#
#		return

		systres = helper.load_obj('systs')
#		for key1 in systres :
#			print key1
#			for key2 in systres[key1] :
#				print key2
#				for key3 in systres[key1][key2] :
#					print key3
		self.make_datacard(systres, 'mm', '++')



#		for syst in systres :
#			for chan in systres[syst] :
#				systres[syst][chan].set_ttwzPredictions()
#				systres[syst][chan].chan_str = chan
#		self.make_datacard(systres, 'em')
#
#		# saving results to file (just for development
#		helper.save_obj(res_syst, 'bla')

#		res_safe = helper.load_obj('bla')
#		self.print_results(res_safe['Normal'])
#
#		self.make_datacard(res_safe, 'em')
#
#		raw_input('ok? ')

		## cout << "=== Going to call makeRatioControlPlots and fillRatios methods..." << endl;
		## if(readHistos(fOutputFileName) != 0) return;
		## 
		## #/* Here I calculate the EWK MC scale factors, but I don't use 
		## #   them in producing any of the fakerate control plots   */
		## bool saveRatioControlPlots = true;
		## makeRatioControlPlots(1, true, saveRatioControlPlots); // El
		## makeRatioControlPlots(2, true, saveRatioControlPlots); // Mu17
		## makeRatioControlPlots(3, true, saveRatioControlPlots); // Mu24_eta2p1
		## 
		## #  /* Here I produce the control plots without re-calculating the
		## #     EWK MC scale factors, but instead using the numbers obtained above. 
		## #     NB: without the previous calls, a SF=1 would be used (FIXME: this is bad design)*/
		## makeRatioControlPlots(1, false, saveRatioControlPlots); // El
		## makeRatioControlPlots(2, false, saveRatioControlPlots); // Mu17
		## makeRatioControlPlots(3, false, saveRatioControlPlots); // Mu24_eta2p1
		## 
		## 
		## bool saveRatioPlots = false;
		## fillRatios(fMuTotData,    fEGData,    0, false, saveRatioPlots);  //make FR plots without applying EWK subtraction
		## fillRatios(fMuTotData,    fEGData,    0, true,  saveRatioPlots);  //make the same applying the EWK subtraction
		## 
		## 
		## fillRatios(fMCBGMuEnr, fMCBGEMEnr, 1, false, false);
		## storeWeightedPred(gRegion[gBaseRegion]);
		## ttG_SR0 = setTTGammaPred(gRegion["SR00"]);
		## cout << "...done ====" << endl;


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

		print '%s SF = data / (WJets + DYJets) = %f / %f = %f' % (chan_str, n_data, n_mc, n_data/n_mc)

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


	def storeWeightedPred(self) :
		foo = 0
		# only nsst, nssl, nzt, nzl for each sample and channel is needed in calculateRatio -> just read those from file and don't loop over tree

		## if(chan != ElMu){
		##         int muelswitch = -1;
		##         if      (chan == Muon) muelswitch = 0;
		## 	else if (chan == Elec) muelswitch = 1;
		## 	S->numbers[reg][chan].nsst = S->tlratios[muelswitch].fntight->GetEntries();
		## 	S->numbers[reg][chan].nssl = S->tlratios[muelswitch].fnloose->GetEntries();
		## 	S->numbers[reg][chan].nzt  = S->tlratios[muelswitch].pntight->GetEntries();
		## 	S->numbers[reg][chan].nzl  = S->tlratios[muelswitch].pnloose->GetEntries();
		## }


	def make_IntPredictions(self, sel) :
		'''oberservation and prediction for different selections'''
		## for now only with one selection

		##################
		# INIT VARIABLES #
		##################

		FR = ROOT.FakeRatios()

		yields = prediction.prediction()

		res = {}
		for ch_str, charge in self.charges.iteritems() :
			res[ch_str] = {}
			res[ch_str]['al'] = result.result('al', charge, 'al')
			res[ch_str]['mm'] = result.result('mm', charge, 'mm')
			res[ch_str]['em'] = result.result('em', charge, 'em')
			res[ch_str]['ee'] = result.result('ee', charge, 'ee')

#		# tight-tight, tight-loose, loose-tight and loose-loose data yields
#		nt2_mm = 0.; nt10_mm = 0.; nt01_mm = 0.; nt0_mm = 0.;
#		nt2_em = 0.; nt10_em = 0.; nt01_em = 0.; nt0_em = 0.;
#		nt2_ee = 0.; nt10_ee = 0.; nt01_ee = 0.; nt0_ee = 0.;

#		# FR Predictions from event-by-event weights (pre stored)
#		npp_mm = 0.; npf_mm = 0.; nfp_mm = 0.; nff_mm = 0.;
#		npp_em = 0.; npf_em = 0.; nfp_em = 0.; nff_em = 0.;
#		npp_ee = 0.; npf_ee = 0.; nfp_ee = 0.; nff_ee = 0.;

#		# OS yields (tight-tight)
#		nt2_ee_BB_os = 0.; nt2_ee_EE_os = 0.; nt2_ee_EB_os = 0.;
#		nt2_em_BB_os = 0.; nt2_em_EE_os = 0.;

		# charge factor: only takes half of the ChMisID prediction if a charge selection is applied
		chargeFactor = 1.
		if sel.charge != 0 : chargeFactor = 0.5

		# rare SM yields (tight-tight)
		nt2_rare_mc_mm = 0.;    nt2_rare_mc_em = 0.;    nt2_rare_mc_ee = 0.;
		nt2_rare_mc_mm_e2 = 0.; nt2_rare_mc_em_e2 = 0.; nt2_rare_mc_ee_e2 = 0.;

		# WZ yields (tight-tight)
		nt2_wz_mc_mm = 0.;    nt2_wz_mc_em = 0.;    nt2_wz_mc_ee = 0.;
		nt2_wz_mc_mm_e2 = 0.; nt2_wz_mc_em_e2 = 0.; nt2_wz_mc_ee_e2 = 0.;

		# all rares, ttW and ttZ yields (tight-tight)
		rares = {}
		rares_npass = {}
#
#		for s in self.samples :
#			if self.samples[s].datamc == 0 : continue
#			yields.rares_mm[s] = 0.; yields.rares_mm_npass[s] = 0;
#			yields.rares_em[s] = 0.; yields.rares_em_npass[s] = 0;
#			yields.rares_ee[s] = 0.; yields.rares_ee_npass[s] = 0;

		nt2_all_mm = {}
		nt2_all_em = {}
		nt2_all_ee = {}
		h_nt2 = {}

		rares_sum = 0.

##		# check if this is faster
##		for s in self.samples :
##			if self.samples[s].datamc == 0 : continue
##			else : scale = self.lumi / self.samples[s].getLumi()
###			print self.sigtree.GetEntries(sel.get_selectionString(False) + ' && SName == ' + s)
###			print s, scale * self.sigtree.GetEntries('PUWeight * HLTSF * ('+sel.get_selectionString(False)+'&& SName == \"'+s+'\")')
##
##			if s == 'WWTo2L2Nu' : continue
##			h_nt2[s] = ROOT.TH1D(s+'_flavor', s + ' flavor', 3, 0., 3.)
##			self.sigtree.Draw('Flavor>>'+s, '%s * PUWeight * HLTSF * (' % (scale) + sel.get_selectionString() + '&& SType == 15 && SName == \"'+s+'\")', 'goff')
###			print 'Flavor>>'+s, '%s * PUWeight * HLTSF * ('%(scale)+sel.get_selectionString(False)+'&& SName == \"'+s+'\")'
##			rares_sum += h_nt2[s].Integral()
##			print s, h_nt2[s].Integral()
###			if event.SType is 15 and event.TLCat == 0 and event.Flavor < 3 :
###				if event.SName == 'WWTo2L2Nu' : continue # TODO: why?
##
##		print '============================='
##		print 'rares sum:', rares_sum

		print sel.get_selectionString()

		#######################
		# ChMisID predictions #
		#######################

		if sel.systflag is 0 :
			print '[status] getting opposite-sign yields..'

			# EM OS
			nt2_em_BB_os = chargeFactor * self.sigtree.GetEntries(sel.get_selectionString((4,0)))
			nt2_em_EE_os = chargeFactor * self.sigtree.GetEntries(sel.get_selectionString((4,1)))

			# EE OS
			nt2_ee_BB_os = chargeFactor * self.sigtree.GetEntries(sel.get_selectionString((5,0)))
			nt2_ee_EB_os = chargeFactor * self.sigtree.GetEntries(sel.get_selectionString((5,4)))
			nt2_ee_EE_os = chargeFactor * self.sigtree.GetEntries(sel.get_selectionString((5,3)))

			(fbb, fbbE) = self.calculateChMisIdProb(self.get_samples('DoubleEle'), 'BB', self.chmid_sf)
			(feb, febE) = self.calculateChMisIdProb(self.get_samples('DoubleEle'), 'EB', self.chmid_sf)
			(fee, feeE) = self.calculateChMisIdProb(self.get_samples('DoubleEle'), 'EE', self.chmid_sf)

	#		(fbb_mc, fbbE_mc) = calculateChMisIdProb(fMCBG, 'BB')
	#		(feb_mc, febE_mc) = calculateChMisIdProb(fMCBG, 'EB')
	#		(fee_mc, feeE_mc) = calculateChMisIdProb(fMCBG, 'EE')

			# Simple error propagation assuming error on number of events is sqrt(N)
			nt2_ee_chmid    = 2*fbb* nt2_ee_BB_os                           + 2*fee*nt2_ee_EE_os                      + 2*feb*nt2_ee_EB_os;
			nt2_ee_chmid_e1 = math.sqrt( 4*fbb*fbb * FR.getEStat2(nt2_ee_BB_os)  + 4*fee*fee * FR.getEStat2(nt2_ee_EE_os)    + 4*feb*feb * FR.getEStat2(nt2_ee_EB_os) )  # stat only
			nt2_ee_chmid_e2 = math.sqrt( 4*fbbE*fbbE * nt2_ee_BB_os*nt2_ee_BB_os + 4*feeE*feeE * nt2_ee_EE_os*nt2_ee_EE_os   + 4*febE*febE * nt2_ee_EB_os*nt2_ee_EB_os + self.ChMisESyst2 * nt2_ee_chmid*nt2_ee_chmid )  # syst only

			nt2_em_chmid    = fbb * nt2_em_BB_os + fee * nt2_em_EE_os;
			nt2_em_chmid_e1 = math.sqrt( fbb*fbb * FR.getEStat2(nt2_em_BB_os)  + fee*fee*FR.getEStat2(nt2_em_EE_os) )
			nt2_em_chmid_e2 = math.sqrt( fbbE*fbbE * nt2_em_BB_os*nt2_em_BB_os + feeE*feeE * nt2_em_EE_os*nt2_em_EE_os + self.ChMisESyst2*nt2_em_chmid*nt2_em_chmid )

			for ch_str, charge in self.charges.iteritems() :
				factor = 1.
				if charge != 0: factor = 0.5
				# TODO: take care about the correct place to apply the chargeFactor to the uncertainties!!
				res[ch_str]['al'].cmid     = factor * nt2_ee_chmid + nt2_em_chmid
				res[ch_str]['em'].cmid     = factor * nt2_em_chmid
				res[ch_str]['ee'].cmid     = factor * nt2_ee_chmid
				res[ch_str]['al'].cmid_err = factor * math.sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2 + nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2)
				res[ch_str]['em'].cmid_err = factor * math.sqrt(nt2_em_chmid_e1*nt2_em_chmid_e1 + nt2_em_chmid_e2*nt2_em_chmid_e2)
				res[ch_str]['ee'].cmid_err = factor * math.sqrt(nt2_ee_chmid_e1*nt2_ee_chmid_e1 + nt2_ee_chmid_e2*nt2_ee_chmid_e2)

		else :
			print '[status] skipped charge mis-ID prediction for systematics'

		################
		# Results tree #
		################

#		results_file = ROOT.TFile(self.path + '/SSDLResults.root', 'RECREATE')
#		results_tree = ROOT.TTree('Results', 'ResultsTree')
#
#		# create branches
#		results_tree.Branch("Weight",      Weight     , "Weight/F"   );
#		results_tree.Branch("SystFlag",    SystFlag   , "SystFlag/I" );
#		results_tree.Branch("SName",       SName);
#		results_tree.Branch("SType",       SType      , "SType/I"    );
#		results_tree.Branch("Run",         Run        , "Run/I"      );
#		results_tree.Branch("LS",          LS         , "LS/I"       );
#		results_tree.Branch("Event",       Event      , "Event/I"    );
#		results_tree.Branch("Flavor",      Flavor     , "Flavor/I"   );
#		results_tree.Branch("Charge",      Charge     , "Charge/I"   );
#		results_tree.Branch("TLCat",       TLCat      , "TLCat/I"    );
#		results_tree.Branch("PassZVeto",   ZVeto      , "PassZVeto/I");
#		results_tree.Branch("HT",          HT         , "HT/F"       );
#		results_tree.Branch("MET",         MET        , "MET/F"      );
#		results_tree.Branch("NJ",          NJ         , "NJ/I"       );
#		results_tree.Branch("NbJ",         NbJ        , "NbJ/I"      );
#		results_tree.Branch("NbJmed",      NbJmed     , "NbJmed/I"   );
#		results_tree.Branch("Mll",         Mll        , "Mll/F"      );
#		results_tree.Branch("pT1",         pT1        , "pT1/F"      );
#		results_tree.Branch("pT2",         pT2        , "pT2/F"      );
#		results_tree.Branch("PFIso1",      PFIso1     , "PFIso1/F"   );
#		results_tree.Branch("PFIso2",      PFIso2     , "PFIso2/F"   );
#		results_tree.Branch("D01",         D01        , "D01/F"      );
#		results_tree.Branch("D02",         D02        , "D02/F"      );
#		results_tree.Branch("Rho",         Rho        , "Rho/F"      );
#		results_tree.Branch("MTLep1",      MTLep1     , "MTLep1/F"   );
#		results_tree.Branch("MTLep2",      MTLep2     , "MTLep2/F"   );
#		results_tree.Branch("dRbJl1",      dRbJl1     , "dRbJl1/F"   );
#		results_tree.Branch("dRbJl2",      dRbJl2     , "dRbJl2/F"   );
#		results_tree.Branch("BetaStar1",   BetaStar1  , "BetaStar1/F");
#		results_tree.Branch("BetaStar2",   BetaStar2  , "BetaStar2/F");
#		results_tree.Branch("BetaStar3",   BetaStar3  , "BetaStar3/F");
#		results_tree.Branch("BetaStar4",   BetaStar4  , "BetaStar4/F");
#		results_tree.Branch("BetaStar5",   BetaStar5  , "BetaStar5/F");
#		results_tree.Branch("NVrtx",       NVrtx      , "NVrtx/I"    );
#		results_tree.Branch("M3",          M3         , "M3/F"       );

		####################################
		# Loop over skimmed SigEvents tree #
		####################################

		last_sample = ''
		nevents = {}
		for s in self.samples : nevents[s] = 0

		print '[status] looping over skimmed tree'
		for i, event in enumerate(self.skimtree) :
			nevents[str(event.SName)] += 1
			if last_sample != str(event.SName) :
				print '[status] processing %s..' % (event.SName)
				last_sample = str(event.SName)

			if not sel.passes_selection(event = event, ttLeptons = False) : continue
			chan = self.get_channelString(int(event.Flavor))

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

#				# EM OS
#				if event.Flavor is 4 :
#					if event.TLCat is 0 : nt2_em_BB_os += chargeFactor
#					if event.TLCat is 1 : nt2_em_EE_os += chargeFactor
#
#				# EE OS
#				if event.Flavor is 5 :
#					if event.TLCat is 0                     : nt2_ee_BB_os += chargeFactor
#					if event.TLCat is 1 or event.TLCat is 2 : nt2_ee_EB_os += chargeFactor
#					if event.TLCat is 3                     : nt2_ee_EE_os += chargeFactor

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

		# END LOOP OVER TREE

		# print number of events per sample in skimmed SigEvents tree
#		for name, n in nevents.iteritems() :
#			print '%12s: %8d' % (name, n)
#		for name, value in rares_mm.iteritems() :

		################
		# Observations #
		################

		for ch_str in res :
			for chan in res[ch_str] :
				res[ch_str][chan].set_observations()

		#####################
		# Fakes predictions #
		#####################

		FR.setNToyMCs(100)
		FR.setAddESyst(0.5)

		# numbers from SSDLPlotter.cc just to check
		mufratio_data = 0.040942; mufratio_data_e = 0.002156;
		mupratio_data = 0.804292; mupratio_data_e = 0.001193;
		elfratio_data = 0.069959; elfratio_data_e = 0.001419;
		elpratio_data = 0.750609; elpratio_data_e = 0.001473;
		FR.setMFRatio(mufratio_data, mufratio_data_e) # set error to pure statistical of ratio
		FR.setEFRatio(elfratio_data, elfratio_data_e)
		FR.setMPRatio(mupratio_data, mupratio_data_e)
		FR.setEPRatio(elpratio_data, elpratio_data_e)

#		# TODO: use these (ratios with ewk subtraction)
#		FR.setMFRatio(self.fpr.MufRatio, self.fpr.MufRatioE) # set error to pure statistical of ratio
#		FR.setEFRatio(self.fpr.ElfRatio, self.fpr.ElfRatioE)
#		FR.setMPRatio(self.fpr.MupRatio, self.fpr.MupRatioE)
#		FR.setEPRatio(self.fpr.ElpRatio, self.fpr.ElpRatioE)

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

		for ch_str in self.charges :

			for chan in res[ch_str] :

				wz   = 0.; wz_staterr2   = 0.;
				ttw  = 0.; ttw_staterr2  = 0.;
				ttz  = 0.; ttz_staterr2  = 0.;
				rare = 0.; rare_staterr2 = 0.;

				for s in rares :
					scale = self.lumi / self.samples[s].getLumi()
					staterr = scale * self.samples[s].getError(rares_npass[s][ch_str][chan])

					print 'sample', s

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

		scale_syst_up = 1.027
		scale_syst_dn = 0.980
		tmass_syst_up = 1.019
		tmass_syst_dn = 0.983
		ltrig_syst    = 1.050
		pdf_syst      = 1.015
		gen_syst      = 1.050

		timestamp = time.asctime()
		print '#========================================================================================='
		print '# Systematics table for ttW analysis, same-sign channel, subchannels'
		print '# Generated on: %s' % (timestamp)
		print '# Copy between the dashed lines for datacard'
		print '#-----------------------------------------------------------------------------------------'
		print 'imax 1'
		print 'jmax 5'
		print 'kmax *'
		print ''
		print 'bin\t\t%s' % (results['Normal'][charge][chan].chan_str)
		##	if (gFullDataBlind)
		##		fOUTSTREAM << Form("observation\t%d\t%d\t%d\t%d\t%d\t%d", 999, 999, 999, 999, 999, 999) << endl;
		##	else
		print 'observation\t%d' % (results['Normal'][charge][chan].obs)
		print '\n'
		print 'bin\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s' % (
			results['Normal'][charge][chan].chan_str,
			results['Normal'][charge][chan].chan_str,
			results['Normal'][charge][chan].chan_str,
			results['Normal'][charge][chan].chan_str,
			results['Normal'][charge][chan].chan_str,
			results['Normal'][charge][chan].chan_str)
		print 'process\t\tttW\t\tttZ\t\tfake\t\tcmid\t\twz\t\trare'
		print 'process\t\t0\t\t1\t\t2\t\t3\t\t4\t\t5'
		print 'rate\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f' % (
			results['Normal'][charge][chan].ttw,
			results['Normal'][charge][chan].ttz,
			results['Normal'][charge][chan].fake,
			results['Normal'][charge][chan].cmid,
			results['Normal'][charge][chan].wz,
			results['Normal'][charge][chan].rare)
		print '\n'
		print '#syst'
		print 'bgUncttz lnN\t-\t\t%5.3f\t\t-\t\t-\t\t-\t\t-' % (1. + results['Normal'][charge][chan].ttz_err  / results['Normal'][charge][chan].ttz )
		print 'bgUncfak lnN\t-\t\t-\t\t%5.3f\t\t-\t\t-\t\t-' % (1. + results['Normal'][charge][chan].fake_err / results['Normal'][charge][chan].fake)
		if chan != 'mm' :
			print 'bgUnccmi lnN\t-\t\t-\t\t-\t\t%5.3f\t\t-\t\t-' % (1. + results['Normal'][charge][chan].cmid_err / results['Normal'][charge][chan].cmid)
		print 'bgUncwz  lnN\t-\t\t-\t\t-\t\t-\t\t%5.3f\t\t-' % (1. + results['Normal'][charge][chan].wz_err   / results['Normal'][charge][chan].wz  )
		print 'bgUncrar lnN\t-\t\t-\t\t-\t\t-\t\t-\t\t%5.3f' % (1. + results['Normal'][charge][chan].rare_err / results['Normal'][charge][chan].rare)
		if 'LepUp' in results and 'LepDown' in results :
			print 'lept     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f' % (
				results['LepDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['LepUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['LepDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['LepUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['LepDown'][charge][chan].wz   / results['Normal'][charge][chan].wz  ,
				results['LepUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz  ,
				results['LepDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
				results['LepUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare)
		else :
			print '[WARNING] LepUp/Down systematic not found!'
		if 'MuUp' in results and 'MuDown' in results :
			print 'muon     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f' % (
				results['MuDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['MuUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['MuDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['MuUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['MuDown'][charge][chan].wz   / results['Normal'][charge][chan].wz  ,
				results['MuUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz  ,
				results['MuDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
				results['MuUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare)
		else :
			print '[WARNING] MuUp/Down systematic not found!'
		if 'ElUp' in results and 'ElDown' in results :
			print 'elec     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f' % (
				results['ElDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['ElUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['ElDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['ElUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['ElDown'][charge][chan].wz   / results['Normal'][charge][chan].wz  ,
				results['ElUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz  ,
				results['ElDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
				results['ElUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare)
		else :
			print '[WARNING] ElUp/Down systematic not found!'
		print 'leptrig  lnN\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f' % (
			ltrig_syst,
			ltrig_syst,
			ltrig_syst,
			ltrig_syst)
		if 'BUp' in results and 'BDown' in results :
			print 'btag     lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f' % (
				results['BDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['BUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['BDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['BUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['BDown'][charge][chan].wz   / results['Normal'][charge][chan].wz,
				results['BUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz,
				results['BDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
				results['BUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare)
		else :
			print '[WARNING] BUp/Down systematic not found!'
		if 'JetUp' in results and 'JetDown' in results :
			print 'jes      lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f' % (
				results['JetDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['JetUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['JetDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['JetUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['JetDown'][charge][chan].wz   / results['Normal'][charge][chan].wz,
				results['JetUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz,
				results['JetDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
				results['JetUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare)
		else :
			print '[WARNING] JetUp/Down systematic not found!'
		##//	fOUTSTREAM << Form("jer      lnN\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t%5.3f\t\t-\t\t-\t\t%5.3f\t\t%5.3f",
		##//					   1.0+(ttwzpreds_plpl["JetSmear"].ttwz_mm-ttwzpreds_plpl["Normal"].ttwz_mm)/ttwzpreds_plpl["Normal"].ttwz_mm,
		##//					   1.0+(ttwzpreds_plpl["JetSmear"].ttwz_mm-ttwzpreds_plpl["Normal"].ttwz_mm)/ttwzpreds_plpl["Normal"].ttwz_mm,
		##//					   1.0+(ttwzpreds_plpl["JetSmear"].wz_mm  -ttwzpreds_plpl["Normal"].wz_mm  )/ttwzpreds_plpl["Normal"].wz_mm,
		##//					   1.0+(ttwzpreds_plpl["JetSmear"].rare_mm-ttwzpreds_plpl["Normal"].rare_mm)/ttwzpreds_plpl["Normal"].rare_mm,
		##//					   1.0+(ttwzpreds_plpl["JetSmear"].ttwz_em-ttwzpreds_plpl["Normal"].ttwz_em)/ttwzpreds_plpl["Normal"].ttwz_em,
		##//					   1.0+(ttwzpreds_plpl["JetSmear"].ttwz_em-ttwzpreds_plpl["Normal"].ttwz_em)/ttwzpreds_plpl["Normal"].ttwz_em,
		##//					   1.0+(ttwzpreds_plpl["JetSmear"].wz_em  -ttwzpreds_plpl["Normal"].wz_em  )/ttwzpreds_plpl["Normal"].wz_em,
		##//					   1.0+(ttwzpreds_plpl["JetSmear"].rare_em-ttwzpreds_plpl["Normal"].rare_em)/ttwzpreds_plpl["Normal"].rare_em,
		##//					   1.0+(ttwzpreds_plpl["JetSmear"].ttwz_ee-ttwzpreds_plpl["Normal"].ttwz_ee)/ttwzpreds_plpl["Normal"].ttwz_ee,
		##//					   1.0+(ttwzpreds_plpl["JetSmear"].ttwz_ee-ttwzpreds_plpl["Normal"].ttwz_ee)/ttwzpreds_plpl["Normal"].ttwz_ee,
		##//					   1.0+(ttwzpreds_plpl["JetSmear"].wz_ee  -ttwzpreds_plpl["Normal"].wz_ee  )/ttwzpreds_plpl["Normal"].wz_ee,
		##//					   1.0+(ttwzpreds_plpl["JetSmear"].rare_ee-ttwzpreds_plpl["Normal"].rare_ee)/ttwzpreds_plpl["Normal"].rare_ee,
		##//					   1.0+(ttwzpreds_mimi["JetSmear"].ttwz_mm-ttwzpreds_mimi["Normal"].ttwz_mm)/ttwzpreds_mimi["Normal"].ttwz_mm,
		##//					   1.0+(ttwzpreds_mimi["JetSmear"].ttwz_mm-ttwzpreds_mimi["Normal"].ttwz_mm)/ttwzpreds_mimi["Normal"].ttwz_mm,
		##//					   1.0+(ttwzpreds_mimi["JetSmear"].wz_mm  -ttwzpreds_mimi["Normal"].wz_mm  )/ttwzpreds_mimi["Normal"].wz_mm,
		##//					   1.0+(ttwzpreds_mimi["JetSmear"].rare_mm-ttwzpreds_mimi["Normal"].rare_mm)/ttwzpreds_mimi["Normal"].rare_mm,
		##//					   1.0+(ttwzpreds_mimi["JetSmear"].ttwz_em-ttwzpreds_mimi["Normal"].ttwz_em)/ttwzpreds_mimi["Normal"].ttwz_em,
		##//					   1.0+(ttwzpreds_mimi["JetSmear"].ttwz_em-ttwzpreds_mimi["Normal"].ttwz_em)/ttwzpreds_mimi["Normal"].ttwz_em,
		##//					   1.0+(ttwzpreds_mimi["JetSmear"].wz_em  -ttwzpreds_mimi["Normal"].wz_em  )/ttwzpreds_mimi["Normal"].wz_em,
		##//					   1.0+(ttwzpreds_mimi["JetSmear"].rare_em-ttwzpreds_mimi["Normal"].rare_em)/ttwzpreds_mimi["Normal"].rare_em,
		##//					   1.0+(ttwzpreds_mimi["JetSmear"].ttwz_ee-ttwzpreds_mimi["Normal"].ttwz_ee)/ttwzpreds_mimi["Normal"].ttwz_ee,
		##//					   1.0+(ttwzpreds_mimi["JetSmear"].ttwz_ee-ttwzpreds_mimi["Normal"].ttwz_ee)/ttwzpreds_mimi["Normal"].ttwz_ee,
		##//					   1.0+(ttwzpreds_mimi["JetSmear"].wz_ee  -ttwzpreds_mimi["Normal"].wz_ee  )/ttwzpreds_mimi["Normal"].wz_ee,
		##//					   1.0+(ttwzpreds_mimi["JetSmear"].rare_ee-ttwzpreds_mimi["Normal"].rare_ee)/ttwzpreds_mimi["Normal"].rare_ee) << endl;
		if 'JetSmearUp' in results and 'JetSmearDown' in results :
			print 'jer      lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f' % (
				results['JetSmearDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['JetSmearUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['JetSmearDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['JetSmearUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['JetSmearDown'][charge][chan].wz   / results['Normal'][charge][chan].wz,
				results['JetSmearUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz,
				results['JetSmearDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
				results['JetSmearUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare)
		else :
			print '[WARNING] JetSmearUp/Down systematic not found!'
		if 'PileupUp' in results and 'PileupDown' in results :
			print 'pu       lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f' % (
				results['PileupDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['PileupUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['PileupDown'][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['PileupUp'  ][charge][chan].ttwz / results['Normal'][charge][chan].ttwz,
				results['PileupDown'][charge][chan].wz   / results['Normal'][charge][chan].wz,
				results['PileupUp'  ][charge][chan].wz   / results['Normal'][charge][chan].wz,
				results['PileupDown'][charge][chan].rare / results['Normal'][charge][chan].rare,
				results['PileupUp'  ][charge][chan].rare / results['Normal'][charge][chan].rare)
		else :
			print '[WARNING] PileupUp/Down systematic not found!'
		##//	fOUTSTREAM << Form("matching lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f",
		##//					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn,
		##//					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn,
		##//					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn,
		##//					   match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn, match_syst_up, match_syst_dn) << endl;
		print 'scale    lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f' % (
			scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn, scale_syst_up, scale_syst_dn)
		print 'tmass    lnN\t%5.3f/%5.3f\t%5.3f/%5.3f\t-\t\t-\t\t%5.3f/%5.3f\t%5.3f/%5.3f' % (
			tmass_syst_up, tmass_syst_dn, tmass_syst_up, tmass_syst_dn, tmass_syst_up, tmass_syst_dn, tmass_syst_up, tmass_syst_dn)
		print 'gen      lnN\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-' % (gen_syst)
		print 'pdf      lnN\t%5.3f\t\t-\t\t-\t\t-\t\t-\t\t-' % (pdf_syst)


	def make_DiffPredictions(self, selections) :
		'''plot observed and predicted distributions of various variables'''

		# create subdirectories for selections

		# init histograms
		histos = {}
		for selection in selections :
			histos[selection.name] = []
			histos[selection.name].append(('pt1', ROOT.TH1D(selection.name + '_pt1', 'Leading Lepton p_{T} (GeV)', 28, 20., 300.)))

		##########
		# RATIOS #
		##########

		# fake ratios

		# charge mis ID probability

		# loop over sigtree
		for event in self.sigtree :
			print '[status] looping over sigtree..'

			


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
