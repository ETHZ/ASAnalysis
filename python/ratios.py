#! /usr/bin/python
import ROOT
import sample
import helper


class ratios :


	def __init__(self, path, samples) :
		self.path = path
		self.ssdlfile = ROOT.TFile.Open(path + '/SSDLYields.root', 'READ')
		self.samples = samples

		# lumi norm
		self.lumi = 19466.
		self.lumi_HLTMu17       = 24.9
		self.lumi_HLTMu24Eta2p1 = 116.
		self.lumi_HLTEl17Jet30  = 23.845


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