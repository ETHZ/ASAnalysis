#! /usr/bin/python
#import os, sys, commands, subprocess, math
import ROOT
import sys
import selection as sel
import sample
import helper

class plotter :
	'''the plotter reads sigtree and produces plots'''

	def __init__(self, path) :
		print '[status] initialize plotter..'
		self.path = path
		self.ssdlfile = ROOT.TFile.Open(path + '/SSDLYields.root', 'READ')
		self.sigtree = self.ssdlfile.Get('SigEvents')
		print '[status] loaded SigEventsTree with %d events' % (self.sigtree.GetEntries())

		ROOT.gSystem.Load('./FakeRatios.so')

		# samples
		self.samples = {}

		# selections
		self.selections = {}

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


	def do_analysis(self, cardfile) :
		print '[statur] starting analysis..'

		# samples
		self.samples = self.readDatacard(cardfile)
		self.read_ngen()

		# selections
		presel = sel.selection(name = 'presel', minNjets = 3, minNbjetsL = 1, minNbjetsM = 1)
		#print presel.get_selectionString()
		self.selections[presel.name] = presel

		EWK_SF = {}
		EWK_SF['el']   = self.get_EWK_SF('el')
		EWK_SF['mu17'] = self.get_EWK_SF('mu17')
		EWK_SF['mu24'] = self.get_EWK_SF('mu24')

		applyEwkSubtr = True
		self.fill_ratios(self.get_samples('SingleDoubleMu'), self.get_samples('DoubleEle'), 0, True, EWK_SF)

		print 'MufRatio: %f +/- %f' % (self.MufRatio, self.MufRatioE)
		print 'ElfRatio: %f +/- %f' % (self.ElfRatio, self.ElfRatioE)
		print 'MupRatio: %f +/- %f' % (self.MupRatio, self.MupRatioE)
		print 'ElpRatio: %f +/- %f' % (self.ElpRatio, self.ElpRatioE)

		print "self.get_fRatio('Muon', 35., 0.5, 0):", self.get_fRatio('Muon', 35., 0.5, 0)
		print "self.get_fRatio('Muon', 55., 0.5, 0):", self.get_fRatio('Muon', 55., 0.5, 0)
		print "self.get_pRatio('Muon', 55., 0     ):", self.get_pRatio('Muon', 55., 0     )

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
		self.make_IntPredictions(presel)
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


	def fill_ratios(self, mu_samples, el_samples, datamc, applyEwkSubtr = False, EWK_SF = {}) :
		print '[status] filling fake and prompt ratio histograms..'

		if datamc is 0 :
			(self.h2_MufRatio   , self.h_MufRatio_pt   , self.h_MufRatio_eta   , self.h_MufRatio_nv   , self.MufRatio   , self.MufRatioE   ) = self.calculateRatio(mu_samples, 'MM', 'SigSup', applyEwkSubtr , EWK_SF)
			(self.h2_ElfRatio   , self.h_ElfRatio_pt   , self.h_ElfRatio_eta   , self.h_ElfRatio_nv   , self.ElfRatio   , self.ElfRatioE   ) = self.calculateRatio(el_samples, 'EE', 'SigSup', applyEwkSubtr , EWK_SF)
			(self.h2_MupRatio   , self.h_MupRatio_pt   , self.h_MupRatio_eta   , self.h_MupRatio_nv   , self.MupRatio   , self.MupRatioE   ) = self.calculateRatio(mu_samples, 'MM', 'ZDecay')
			(self.h2_ElpRatio   , self.h_ElpRatio_pt   , self.h_ElpRatio_eta   , self.h_ElpRatio_nv   , self.ElpRatio   , self.ElpRatioE   ) = self.calculateRatio(el_samples, 'EE', 'ZDecay')

		else :
			(self.h2_MufRatio_MC, self.h_MufRatio_pt_MC, self.h_MufRatio_eta_MC, self.h_MufRatio_nv_MC, self.MufRatio_MC, self.MufRatioE_MC) = self.calculateRatio(mu_samples, 'MM', 'SigSup', applyEwkSubtr , EWK_SF)
			(self.h2_ElfRatio_MC, self.h_ElfRatio_pt_MC, self.h_ElfRatio_eta_MC, self.h_ElfRatio_nv_MC, self.ElfRatio_MC, self.ElfRatioE_MC) = self.calculateRatio(el_samples, 'EE', 'SigSup', applyEwkSubtr , EWK_SF)
			(self.h2_MupRatio_MC, self.h_MupRatio_pt_MC, self.h_MupRatio_eta_MC, self.h_MupRatio_nv_MC, self.MupRatio_MC, self.MupRatioE_MC) = self.calculateRatio(mu_samples, 'MM', 'ZDecay')
			(self.h2_ElpRatio_MC, self.h_ElpRatio_pt_MC, self.h_ElpRatio_eta_MC, self.h_ElpRatio_nv_MC, self.ElpRatio_MC, self.ElpRatioE_MC) = self.calculateRatio(el_samples, 'EE', 'ZDecay')


	def calculateRatio(self, samples, chan_str, fp, applyEwkSubtr = False, EWK_SF = {}) :
		# - sets up ratio histos
		# - calls getPassedTotal / get_TightLoose to get ntight and nloose histos

		samples_ewk = []
		samples_ewk.append('WJets')
		samples_ewk.append('DYJets')

		(h2_ntight, h2_nloose, h_ntight_nv, h_nloose_nv) = self.get_TightLoose(samples, chan_str, fp)

		if fp is 'SigSup' and applyEwkSubtr :

			(h2_ntight_ewk, h2_nloose_ewk, h_ntight_nv_ewk, h_nloose_nv_ewk) = self.get_TightLoose(samples_ewk, chan_str, fp)

			if chan_str is 'EE' :
				scale = EWK_SF['el'] * self.lumi_HLTEl17Jet30 / self.lumi  # MC samples are scaled to self.lumi. Rescale them now to prescaled triggers and apply correction factor
				h2_ntight  .Add(h2_ntight_ewk  , (-1.) * scale)
				h2_nloose  .Add(h2_nloose_ewk  , (-1.) * scale)
				h_ntight_nv.Add(h_ntight_nv_ewk, (-1.) * scale)
				h_nloose_nv.Add(h_nloose_nv_ewk, (-1.) * scale)

			if chan_str is 'MM' :
				for pt_bin in range(1, h2_ntight.GetXaxis().GetNbins()+1) :
					for eta_bin in range(1, h2_ntight.GetYaxis().GetNbins()+1) :
						bin = h2_ntight.GetBin(pt_bin, eta_bin)
						# MC samples are scaled to self.lumi. Rescale them now to prescaled triggers and apply correction factor
						if (h2_ntight.GetXaxis().GetBinCenter(pt_bin) > 25.) and (h2_ntight.GetYaxis().GetBinCenter(eta_bin) < 2.1) :
							scale = EWK_SF['mu24'] * self.lumi_HLTMu24Eta2p1 / self.lumi
						else :
							scale = EWK_SF['mu17'] * self.lumi_HLTMu17 / self.lumi
#						print 'pt bin (%d) center: %f.3, eta bin (%d) center: %f.3, scale: %f' % (pt_bin, h2_ntight.GetXaxis().GetBinCenter(pt_bin), eta_bin, h2_ntight.GetYaxis().GetBinCenter(eta_bin), scale)
						h2_ntight.AddBinContent(bin, (-1.) * scale * h2_ntight_ewk.GetBinContent(pt_bin, eta_bin))
						h2_nloose.AddBinContent(bin, (-1.) * scale * h2_nloose_ewk.GetBinContent(pt_bin, eta_bin))

		h2_ratio   = h2_ntight  .Clone('h2_'+chan_str+'_'+fp+'_ratio_vs_pt_eta')
		h_ratio_nv = h_ntight_nv.Clone('h_'+chan_str+'_'+fp+'_ratio_vs_nv')
		h2_ratio  .Divide(h2_ntight  , h2_nloose  , 1., 1., 'B')
		h_ratio_nv.Divide(h_ntight_nv, h_nloose_nv, 1., 1., 'B')

		h_nloose_pt  = h2_nloose.ProjectionX()
		h_nloose_eta = h2_nloose.ProjectionY()
		h_ntight_pt  = h2_ntight.ProjectionX()
		h_ntight_eta = h2_ntight.ProjectionY()

		h_ratio_pt  = h_ntight_pt .Clone('h_'+chan_str+'_'+fp+'_ratio_vs_pt')
		h_ratio_eta = h_ntight_eta.Clone('h_'+chan_str+'_'+fp+'_ratio_vs_eta')
		h_ratio_pt .Divide(h_ntight_pt , h_nloose_pt , 1., 1., 'B')
		h_ratio_eta.Divide(h_ntight_eta, h_nloose_eta, 1., 1., 'B')

		# flat ratios
		ntight = h2_ntight.Integral()
		nloose = h2_nloose.Integral()
		(ratio, ratioe) = helper.ratioWithBinomErrors(ntight, nloose)

		return (h2_ratio, h_ratio_pt, h_ratio_eta, h_ratio_nv, ratio, ratioe)
		## 	TString name = "";
		## 	if (fp == SigSup) name += "F";
		## 	if (fp == ZDecay) name += "P";
		## 	name += h_2d->GetName();
		## 
		## 	/*
		## 	for(size_t i = 0; i < samples.size(); ++i){
		## 		int sample = samples[i];
		## 		name += "_";
		## 		name += fSamples[sample]->sname;
		## 	}
		## 	*/
		## 	//FIXME: add something to distinguish between data and mc plots
		## 	//if() name += "_data";
		## 
		## 	if (applyEwkSubtr) name += "_EWKCorrected";
		## 	if(output){
		## //	if (fp == SigSup) {
		## 		fOutputSubDir = "Ratios/";
		## 		printObject(h_2d,     TString("Ratio")    + name, "colz text");
		## 		printObject(h_pt,     TString("RatioPt")  + name, "PE1");
		## 		printObject(h_eta,    TString("RatioEta") + name, "PE1");
		## 		printObject(H_ntight, TString("NTight")   + name, "colz text");
		## 		printObject(H_nloose, TString("NLoose")   + name, "colz text");
		## 		fOutputSubDir = "";
		## 	}


	def get_TightLoose(self, samples, chan_str, fp) :
		## chan_str: 'MM', 'EE'

		for i, s in enumerate(samples) :
			scale = 1.
			if self.samples[s].datamc != 0 : scale = self.lumi / self.samples[s].getLumi()
			if fp is 'SigSup' :
				h2_ntight   = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_fNTight').Clone()
				h2_nloose   = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_fNLoose').Clone()
				h_ntight_nv = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_fNTight_nv').Clone()
				h_nloose_nv = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_fNLoose_nv').Clone()
			elif fp is 'ZDecay' :
				h2_ntight   = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_pNTight')   .Clone()
				h2_nloose   = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_pNLoose')   .Clone()
				h_ntight_nv = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_pNTight_nv').Clone()
				h_nloose_nv = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_pNLoose_nv').Clone()
			if i is 0 :
				h2_passed   = h2_ntight
				h2_total    = h2_nloose
				h_passed_nv = h_ntight_nv
				h_total_nv  = h_nloose_nv
				h2_passed  .Scale(scale)
				h2_total   .Scale(scale)
				h_passed_nv.Scale(scale)
				h_total_nv .Scale(scale)
			else :
				h2_passed  .Add(h2_ntight  , scale)
				h2_total   .Add(h2_nloose  , scale)
				h_passed_nv.Add(h_ntight_nv, scale)
				h_total_nv .Add(h_nloose_nv, scale)

		return (h2_passed, h2_total, h_passed_nv, h_total_nv)
		## 	TString name = "";
		## 	for(size_t i = 0; i < samples.size(); ++i){
		## 		int sample = samples[i];
		## 		if(i > 0) name += "_";
		## 		name += fSamples[sample]->sname;
		## 	}
		## 	/*
		## 	if(output){
		## 		printObject(h_passed, TString("Passed") + name, "colz");
		## 		printObject(h_total,  TString("Total")  + name, "colz");
		## 	}
		## 	*/	
		## }


	def get_fRatio(self, chan, pt, eta, datamc) :
		mu_flatout = 40.
		el_flatout = 40.
		eta = abs(eta)

		if chan is 'Muon' :
			if datamc is 0 : histo = self.h2_MufRatio
			else           : histo = self.h2_MufRatio_MC
			flatout = mu_flatout
		elif chan is 'Elec' :
			if datamc is 0 : histo = self.h2_ElfRatio
			else           : histo = self.h2_ElfRatio_MC
			flatout = el_flatout

##		if histo.GetEntries() is 0   ## TODO: check if histo exists

		if chan != 'Muon' and chan != 'Elec' : print chan

		if pt >= flatout :
			bin = histo.FindBin(flatout-1., eta)
			return histo.GetBinContent(bin)

		bin = histo.FindBin(pt, eta)
		return histo.GetBinContent(bin)


	def get_pRatio(self, chan, pt, datamc) :
		if chan is 'Muon' :
			if datamc is 0 : histo = self.h_MupRatio_pt
			else           : histo = self.h_MupRatio_pt_MC
		elif chan is 'Elec' :
			if datamc is 0 : histo = self.h_ElpRatio_pt
			else           : histo = self.h_ElpRatio_pt_MC

##		if histo.GetEntries() is 0   ## TODO: check if histo exists

		bin = histo.FindBin(pt)
		return histo.GetBinContent(bin)


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

		foo = 0

		# this gets the flat ratio numbers (w/o ewk subtraction)

		##	///////////////////////////////////////////////////////////////////////////////////
		##	// RATIOS /////////////////////////////////////////////////////////////////////////
		##	///////////////////////////////////////////////////////////////////////////////////
		##	float mufratio_data(0.),  mufratio_data_e(0.);
		##	float mupratio_data(0.),  mupratio_data_e(0.);
		##	float elfratio_data(0.),  elfratio_data_e(0.);
		##	float elpratio_data(0.),  elpratio_data_e(0.);
		##
		##	calculateRatio(fMuData, Muon, SigSup, mufratio_data, mufratio_data_e);
		##	calculateRatio(fMuData, Muon, ZDecay, mupratio_data, mupratio_data_e);
		##
		##	calculateRatio(fEGData, Elec, SigSup, elfratio_data, elfratio_data_e);
		##	calculateRatio(fEGData, Elec, ZDecay, elpratio_data, elpratio_data_e);

		# the numbers are calculated earlier with ewk subtraction and stored in self.MufRatio, ...

		##	FakeRatios *FR = new FakeRatios();
		#ROOT.gSystem.Load('./FakeRatios.so')
		FR = ROOT.FakeRatios()

		# init variables
		nt2_mm = 0.; nt10_mm = 0.; nt01_mm = 0.; nt0_mm = 0.;
		nt2_em = 0.; nt10_em = 0.; nt01_em = 0.; nt0_em = 0.;
		nt2_ee = 0.; nt10_ee = 0.; nt01_ee = 0.; nt0_ee = 0.;

		# FR Predictions from event-by-event weights (pre stored)
		npp_mm = 0.; npf_mm = 0.; nfp_mm = 0.; nff_mm = 0.;
		npp_em = 0.; npf_em = 0.; nfp_em = 0.; nff_em = 0.;
		npp_ee = 0.; npf_ee = 0.; nfp_ee = 0.; nff_ee = 0.;

		# OS yields
		nt2_ee_BB_os = 0.; nt2_ee_EE_os = 0.; nt2_ee_EB_os = 0.;
		nt2_em_BB_os = 0.; nt2_em_EE_os = 0.;

		# TODO
		chargeFactor = 1.
		if sel.charge != 0 : chargeFactor = 0.5

		# rare SM yields
		nt2_rare_mc_mm = 0.;    nt2_rare_mc_em = 0.;    nt2_rare_mc_ee = 0.;
		nt2_rare_mc_mm_e2 = 0.; nt2_rare_mc_em_e2 = 0.; nt2_rare_mc_ee_e2 = 0.;

		nt2_wz_mc_mm = 0.;    nt2_wz_mc_em = 0.;    nt2_wz_mc_ee = 0.;
		nt2_wz_mc_mm_e2 = 0.; nt2_wz_mc_em_e2 = 0.; nt2_wz_mc_ee_e2 = 0.;

		rareMapMM = {}
		rareMapEM = {}
		rareMapEE = {}

		rareMapMM_npass = {}
		rareMapEM_npass = {}
		rareMapEE_npass = {}

		for s in self.samples :
			rareMapMM[s] = 0.; rareMapMM_npass[s] = 0;
			rareMapEM[s] = 0.; rareMapEM_npass[s] = 0;
			rareMapEE[s] = 0.; rareMapEE_npass[s] = 0;

		last_sample = ''

		for i, event in enumerate(self.sigtree) :
			if last_sample != str(event.SName) :
				print '[status] processing %s..' % (event.SName)
				last_sample = str(event.SName)
#			if i%100000 is 0 : print '[status] processing event %d' % (i)

			if event.SType > 2 : continue # only data

			if not sel.passes_selection(event, False, True, True) : continue
			chan = self.get_channelString(int(event.Flavor))

			# GET ALL DATA EVENTS
			if event.SType < 3 :
				if event.Flavor < 3 :
					if sel.charge != 0 and event.Charge != sel.charge : continue
					if chan is 'ElMu' :
						f1 = self.get_fRatio('Muon', event.pT1, event.eta1, 0)
						f2 = self.get_fRatio('Elec', event.pT2, event.eta2, 0)
						p1 = self.get_pRatio('Muon', event.pT1, 0)
						p2 = self.get_pRatio('Elec', event.pT2, 0)
					else :
						f1 = self.get_fRatio(chan, event.pT1, event.eta1, 0)
						f2 = self.get_fRatio(chan, event.pT2, event.eta2, 0)
						p1 = self.get_pRatio(chan, event.pT1, 0)
						p2 = self.get_pRatio(chan, event.pT2, 0)

					npp = FR.getWpp(event.TLCat, f1, f2, p1, p2)
					npf = FR.getWpf(event.TLCat, f1, f2, p1, p2)
					nfp = FR.getWfp(event.TLCat, f1, f2, p1, p2)
					nff = FR.getWff(event.TLCat, f1, f2, p1, p2)

					# MM
					if event.Flavor is 0 :
						npp_mm += npp;
						npf_mm += npf;
						nfp_mm += nfp;
						nff_mm += nff;
						if event.TLCat is 0 : nt2_mm  += 1
						if event.TLCat is 1 : nt10_mm += 1
						if event.TLCat is 2 : nt01_mm += 1
						if event.TLCat is 3 : nt0_mm  += 1

					# EM
					if event.Flavor is 1 :
						npp_em += npp
						npf_em += npf
						nfp_em += nfp
						nff_em += nff
						if event.TLCat is 0 : nt2_em  += 1
					 	if event.TLCat is 1 : nt10_em += 1
					 	if event.TLCat is 2 : nt01_em += 1
					 	if event.TLCat is 3 : nt0_em  += 1

					# EE
					if event.Flavor is 2 :
						npp_ee += npp
						npf_ee += npf
						nfp_ee += nfp
						nff_ee += nff
						if event.TLCat is 0 : nt2_ee  += 1
						if event.TLCat is 1 : nt10_ee += 1
						if event.TLCat is 2 : nt01_ee += 1
						if event.TLCat is 3 : nt0_ee  += 1

				# EM OS
				if event.Flavor is 4 :
					if event.TLCat is 0 : nt2_em_BB_os += chargeFactor
					if event.TLCat is 1 : nt2_em_EE_os += chargeFactor

				# EE OS
				if event.Flavor is 5 :
					if event.TLCat is 0                     : nt2_ee_BB_os += chargeFactor
					if event.TLCat is 1 or event.TLCat is 2 : nt2_ee_EB_os += chargeFactor
					if event.TLCat is 3                     : nt2_ee_EE_os += chargeFactor

			# GET RARE MC EVENTS
			if event.SType is 15 and event.TLCat == 0 and event.Flavor < 3 :
				if event.SName == 'WWTo2L2Nu' : continue # TODO: why?
				if sel.charge != 0 and event.Charge != sel.charge : continue
				scale = event.PUWeight * event.HLTSF * self.lumi / self.samples[str(event.SName)].getLumi()

				if event.Flavor is 0 :
					rareMapMM      [str(event.SName)] += scale
					rareMapMM_npass[str(event.SName)] += 1

				if event.Flavor is 1 :
					rareMapEM      [str(event.SName)] += scale
					rareMapEM_npass[str(event.SName)] += 1

				if event.Flavor is 2 :
					rareMapEE      [str(event.SName)] += scale
					rareMapEE_npass[str(event.SName)] += 1

		# end of sigevents loop

		# PRINTOUT
		print "-------------------------------------------------------------------------------------------------------------"
		print "                 |               Mu/Mu               |                E/Mu               |                E/E                ||"
		print "         YIELDS  |   Ntt  |   Ntl  |   Nlt  |   Nll  |   Ntt  |   Ntl  |   Nlt  |   Nll  |   Ntt  |   Ntl  |   Nlt  |   Nll  ||"
		print "-------------------------------------------------------------------------------------------------------------"
		print "%16s & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f & %6.0f" % ("Data", nt2_mm, nt10_mm, nt01_mm, nt0_mm, nt2_em, nt10_em, nt01_em, nt0_em, nt2_ee, nt10_ee, nt01_ee, nt0_ee)


		# numbers from SSDLPlotter.cc just to check
		mufratio_data = 0.040942; mufratio_data_e = 0.002156;
		mupratio_data = 0.804292; mupratio_data_e = 0.001193;
		elfratio_data = 0.069959; elfratio_data_e = 0.001419;
		elpratio_data = 0.750609; elpratio_data_e = 0.001473;

		FR.setNToyMCs(100)
		FR.setAddESyst(0.5)

		# numbers from SSDLPlotter.cc just to check
		FR.setMFRatio(mufratio_data, mufratio_data_e) # set error to pure statistical of ratio
		FR.setEFRatio(elfratio_data, elfratio_data_e)
		FR.setMPRatio(mupratio_data, mupratio_data_e)
		FR.setEPRatio(elpratio_data, elpratio_data_e)

#		# TODO: use these (ratios with ewk subtraction)
#		FR.setMFRatio(self.MufRatio, self.MufRatioE) # set error to pure statistical of ratio
#		FR.setEFRatio(self.ElfRatio, self.ElfRatioE)
#		FR.setMPRatio(self.MupRatio, self.MupRatioE)
#		FR.setEPRatio(self.ElpRatio, self.ElpRatioE)

		FR.setMMNtl(nt2_mm, nt10_mm, nt01_mm, nt0_mm)
		FR.setEENtl(nt2_ee, nt10_ee, nt01_ee, nt0_ee)
		FR.setEMNtl(nt2_em, nt10_em, nt01_em, nt0_em)

		nF_mm = npf_mm + nfp_mm + nff_mm
		nF_em = npf_em + nfp_em + nff_em
		nF_ee = npf_ee + nfp_ee + nff_ee
		nSF   = npf_mm + nfp_mm + npf_em + nfp_em + npf_ee + nfp_ee
		nDF   = nff_mm + nff_em + nff_ee
		nF    = nF_mm + nF_em + nF_ee

		print "  Fake Predictions:"
		print "------------------------------------------------------------------------------------------"
		print "                 |          Mu/Mu        |         El/El         |          El/Mu        |"
		print "------------------------------------------------------------------------------------------"
		print " Npp             |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
			npp_mm, FR.getMMNppEStat(), self.FakeESyst*npp_mm,
			npp_ee, FR.getEENppEStat(), self.FakeESyst*npp_ee, 
			npp_em, FR.getEMNppEStat(), self.FakeESyst*npp_em)
		print " Npf             |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
			npf_mm, FR.getMMNpfEStat(), self.FakeESyst*npf_mm,
			npf_ee, FR.getEENpfEStat(), self.FakeESyst*npf_ee, 
			npf_em, FR.getEMNpfEStat(), self.FakeESyst*npf_em)
		print " Nfp             |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
			nfp_mm, FR.getMMNfpEStat(), self.FakeESyst*nfp_mm,
			nfp_ee, FR.getEENfpEStat(), self.FakeESyst*nfp_ee, 
			nfp_em, FR.getEMNfpEStat(), self.FakeESyst*nfp_em)
		print " Nff             |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
			nff_mm, FR.getMMNffEStat(), self.FakeESyst*nff_mm,
			nff_ee, FR.getEENffEStat(), self.FakeESyst*nff_ee, 
			nff_em, FR.getEMNffEStat(), self.FakeESyst*nff_em)
		print "------------------------------------------------------------------------------------------"
		print " Total Fakes     |", "%5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f | %5.1f +/- %5.1f +/- %5.1f |" % (
			nF_mm, FR.getMMTotEStat(), self.FakeESyst*nF_mm,
			nF_ee, FR.getEETotEStat(), self.FakeESyst*nF_ee, 
			nF_em, FR.getEMTotEStat(), self.FakeESyst*nF_em)
		print "------------------------------------------------------------------------------------------"
		print " (Value +/- E_stat +/- E_syst) "
		print "//////////////////////////////////////////////////////////////////////////////////////////"


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
