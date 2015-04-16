#! /usr/bin/python
import ROOT
import sample
import helper
import ttvplot
import ttvStyle


class ratios :


	def __init__(self, path, samples) :
		self.path = path
		if not self.path.endswith('/') : self.path += '/'
		self.ssdlfile = ROOT.TFile.Open(self.path + 'SSDLYields.root', 'READ')
		self.samples = samples

		# lumi norm
		self.lumi = 19466.
		self.lumi_HLTMu17       = 24.9
		self.lumi_HLTMu24Eta2p1 = 116.
		self.lumi_HLTEl17Jet30  = 23.845


	def fill_ratios(self, mu_samples, el_samples, datamc, applyEwkSubtr = False) :
		EWK_SF = {}
		if applyEwkSubtr :
			if datamc is 0 :
				mu17_samples = filter(lambda sample : sample.startswith('DoubleMu'), mu_samples)
				mu24_samples = filter(lambda sample : sample.startswith('SingleMu'), mu_samples)
				EWK_SF['el']   = self.get_EWK_SF(  el_samples, 'el'  )
				EWK_SF['mu17'] = self.get_EWK_SF(mu17_samples, 'mu17')
				EWK_SF['mu24'] = self.get_EWK_SF(mu24_samples, 'mu24')

			else :
				EWK_SF['mu_mc'] = self.get_EWK_SF(mu_samples, 'mu', datamc = datamc)
				EWK_SF['el_mc'] = self.get_EWK_SF(el_samples, 'el', datamc = datamc)

		print '[status] filling fake and prompt ratio histograms..'
		print '         datamc: %d' % datamc
		print '         mu samples: %s' % ', '.join(mu_samples)
		print '         el samples: %s' % ', '.join(el_samples)

		if datamc is 0 :
			(self.h2_MufRatio   , self.h_MufRatio_pt   , self.h_MufRatio_eta   , self.h_MufRatio_nv   , self.MufRatio   , self.MufRatioE   ) = self.calculateRatio(mu_samples, 'MM', 'SigSup', applyEwkSubtr , EWK_SF)
			(self.h2_ElfRatio   , self.h_ElfRatio_pt   , self.h_ElfRatio_eta   , self.h_ElfRatio_nv   , self.ElfRatio   , self.ElfRatioE   ) = self.calculateRatio(el_samples, 'EE', 'SigSup', applyEwkSubtr , EWK_SF)
			(self.h2_MupRatio   , self.h_MupRatio_pt   , self.h_MupRatio_eta   , self.h_MupRatio_nv   , self.MupRatio   , self.MupRatioE   ) = self.calculateRatio(mu_samples, 'MM', 'ZDecay')
			(self.h2_ElpRatio   , self.h_ElpRatio_pt   , self.h_ElpRatio_eta   , self.h_ElpRatio_nv   , self.ElpRatio   , self.ElpRatioE   ) = self.calculateRatio(el_samples, 'EE', 'ZDecay')
			print '         MufRatio: %f +/- %f' % (self.MufRatio, self.MufRatioE)
			print '         ElfRatio: %f +/- %f' % (self.ElfRatio, self.ElfRatioE)
			print '         MupRatio: %f +/- %f' % (self.MupRatio, self.MupRatioE)
			print '         ElpRatio: %f +/- %f' % (self.ElpRatio, self.ElpRatioE)

		else :
#			mu_samples_fRatio = mu_samples
#			el_samples_fRatio = el_samples
#			if applyEwkSubtr :
#				mu_samples_fRatio = filter(lambda sample : sample != 'WJets' and 'DYJets' not in sample, mu_samples)
#				el_samples_fRatio = filter(lambda sample : sample != 'WJets' and 'DYJets' not in sample, el_samples)

#			(self.h2_MufRatio_MC, self.h_MufRatio_pt_MC, self.h_MufRatio_eta_MC, self.h_MufRatio_nv_MC, self.MufRatio_MC, self.MufRatioE_MC) = self.calculateRatio(mu_samples_fRatio, 'MM', 'SigSup')
#			(self.h2_ElfRatio_MC, self.h_ElfRatio_pt_MC, self.h_ElfRatio_eta_MC, self.h_ElfRatio_nv_MC, self.ElfRatio_MC, self.ElfRatioE_MC) = self.calculateRatio(el_samples_fRatio, 'EE', 'SigSup')
			(self.h2_MufRatio_MC, self.h_MufRatio_pt_MC, self.h_MufRatio_eta_MC, self.h_MufRatio_nv_MC, self.MufRatio_MC, self.MufRatioE_MC) = self.calculateRatio(mu_samples, 'MM', 'SigSup', applyEwkSubtr , EWK_SF, datamc)
			(self.h2_ElfRatio_MC, self.h_ElfRatio_pt_MC, self.h_ElfRatio_eta_MC, self.h_ElfRatio_nv_MC, self.ElfRatio_MC, self.ElfRatioE_MC) = self.calculateRatio(el_samples, 'EE', 'SigSup', applyEwkSubtr , EWK_SF, datamc)
			(self.h2_MupRatio_MC, self.h_MupRatio_pt_MC, self.h_MupRatio_eta_MC, self.h_MupRatio_nv_MC, self.MupRatio_MC, self.MupRatioE_MC) = self.calculateRatio(mu_samples, 'MM', 'ZDecay', datamc = datamc)
			(self.h2_ElpRatio_MC, self.h_ElpRatio_pt_MC, self.h_ElpRatio_eta_MC, self.h_ElpRatio_nv_MC, self.ElpRatio_MC, self.ElpRatioE_MC) = self.calculateRatio(el_samples, 'EE', 'ZDecay', datamc = datamc)
			print '         MufRatio: %f +/- %f' % (self.MufRatio_MC, self.MufRatioE_MC)
			print '         ElfRatio: %f +/- %f' % (self.ElfRatio_MC, self.ElfRatioE_MC)
			print '         MupRatio: %f +/- %f' % (self.MupRatio_MC, self.MupRatioE_MC)
			print '         ElpRatio: %f +/- %f' % (self.ElpRatio_MC, self.ElpRatioE_MC)


	def calculateRatio(self, samples, chan_str, fp, applyEwkSubtr = False, EWK_SF = {}, datamc = 0) :
		'''
		- sets up ratio histos
		- calls get_TightLoose to get ntight and nloose histos
		'''

		samples_ewk = []
		samples_ewk += sample.sample.get_samples('WJets' , self.samples)
		samples_ewk += sample.sample.get_samples('DYJets', self.samples)

		(h2_ntight, h2_nloose, h_ntight_nv, h_nloose_nv) = self.get_TightLoose(samples, chan_str, fp)

		if fp is 'SigSup' and applyEwkSubtr :

			(h2_ntight_ewk, h2_nloose_ewk, h_ntight_nv_ewk, h_nloose_nv_ewk) = self.get_TightLoose(samples_ewk, chan_str, fp)

			if chan_str is 'EE' and datamc == 0 :
				scale = EWK_SF['el'] * self.lumi_HLTEl17Jet30 / self.lumi  # MC samples are scaled to self.lumi. Rescale them now to prescaled triggers and apply correction factor
				h2_ntight  .Add(h2_ntight_ewk  , (-1.) * scale)
				h2_nloose  .Add(h2_nloose_ewk  , (-1.) * scale)
				h_ntight_nv.Add(h_ntight_nv_ewk, (-1.) * scale)
				h_nloose_nv.Add(h_nloose_nv_ewk, (-1.) * scale)

			if chan_str is 'MM' and datamc == 0 :
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

			if datamc != 0 :
				ewk_str = ''
				if chan_str is 'EE' : ewk_str = 'el_mc'
				if chan_str is 'MM' : ewk_str = 'mu_mc'
				scale = EWK_SF[ewk_str]
				h2_ntight  .Add(h2_ntight_ewk  , (-1.) * scale)
				h2_nloose  .Add(h2_nloose_ewk  , (-1.) * scale)
				h_ntight_nv.Add(h_ntight_nv_ewk, (-1.) * scale)
				h_nloose_nv.Add(h_nloose_nv_ewk, (-1.) * scale)

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


	def get_TightLoose(self, samples, chan_str, fp, sf = 1., sf_dict = {}) :
		## chan_str: 'MM', 'EE'

		for i, s in enumerate(samples) :
			scale = 1.
			if self.samples[s].datamc != 0 :
				sf_tmp = 1.
				if s in sf_dict : sf_tmp = sf_dict[s]
				scale = self.lumi / self.samples[s].getLumi() * sf * sf_tmp
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


	def get_EWK_SF(self, samples_data, chan_str, datamc = 0, do_fit = False) :
		print '[status] calculating EWK scale factor for %s trigger..' % (chan_str)
		samples_wjets = []
		samples_zjets = []
		samples_qcd = [] # TODO
#		if chan_str is 'el'   : samples_data = self.get_samples('DoubleEle')
#		if chan_str is 'mu17' : samples_data = self.get_samples('DoubleMu')
#		if chan_str is 'mu24' : samples_data = self.get_samples('SingleMu')
		samples_wjets = sample.sample.get_samples('WJets' , self.samples)
		samples_zjets = sample.sample.get_samples('DYJets', self.samples)
		samples_qcd   = sample.sample.get_samples('QCD'   , self.samples)
#		samples_qcd = ['MuEnr15','QCDEM30', 'QCDEM80', 'QCDEM250', 'QCDEM350']

		(h_ntight_data , h_nloose_data ) = self.get_fRatioPlots(samples_data , chan_str, 'MT_MET30', datamc)
		(h_ntight_wjets, h_nloose_wjets) = self.get_fRatioPlots(samples_wjets, chan_str, 'MT_MET30', datamc)
		(h_ntight_zjets, h_nloose_zjets) = self.get_fRatioPlots(samples_zjets, chan_str, 'MT_MET30', datamc)
		(h_ntight_qcd  , h_nloose_qcd  ) = self.get_fRatioPlots(samples_qcd  , chan_str, 'MT_MET30', datamc)

		if do_fit :
#			(h_ntight_data , h_nloose_data ) = self.get_fRatioPlots(samples_data , chan_str, 'MT', datamc)
#			(h_ntight_wjets, h_nloose_wjets) = self.get_fRatioPlots(samples_wjets, chan_str, 'MT', datamc)
#			(h_ntight_zjets, h_nloose_zjets) = self.get_fRatioPlots(samples_zjets, chan_str, 'MT', datamc)
#			(h_ntight_qcd  , h_nloose_qcd  ) = self.get_fRatioPlots(samples_qcd  , chan_str, 'MT', datamc)
			h_nloose_data  = self.get_fRatioPlots(samples_data , chan_str, 'LepEta', datamc)[1]
			h_nloose_wjets = self.get_fRatioPlots(samples_wjets, chan_str, 'LepEta', datamc)[1]
			h_nloose_zjets = self.get_fRatioPlots(samples_zjets, chan_str, 'LepEta', datamc)[1]
			h_nloose_qcd   = self.get_fRatioPlots(samples_qcd  , chan_str, 'LepEta', datamc)[1]
			ratios = {}
			ratios['tight'] = self.fit_fRatioPlots(h_ntight_data, h_ntight_wjets, h_ntight_zjets, h_ntight_qcd, individual_fit = True)
			ratios['loose'] = self.fit_fRatioPlots(h_nloose_data, h_nloose_wjets, h_nloose_zjets, h_nloose_qcd, individual_fit = True)
			print ratios
			return ratios

		else :
			bin_min = h_ntight_data.FindBin(60.)
			bin_max = h_ntight_data.FindBin(90.)-1

			n_data = h_ntight_data.Integral(bin_min, bin_max)
			n_mc   = h_ntight_wjets.Integral(bin_min, bin_max) + h_ntight_zjets.Integral(bin_min, bin_max)

			ratio = n_data / n_mc

			print '         SF = data / (WJets + DYJets) = %8.1f / %8.1f = %4.2f' % (n_data, n_mc, ratio)

			return ratio


	def get_fRatioPlots(self, samples, chan_str, ratiovar, datamc = 0) :
		'''gets ntight and loose histograms for various variables'''
		## chan_str = 'mu', 'el', 'mu17', 'mu24'
		lumi = self.lumi
		if datamc == 0 :
			if chan_str == 'el'   : lumi = self.lumi_HLTEl17Jet30
			if chan_str == 'mu17' : lumi = self.lumi_HLTMu17
			if chan_str == 'mu24' : lumi = self.lumi_HLTMu24Eta2p1
		for i, s in enumerate(samples) :
			scale = lumi / self.samples[s].getLumi()
			if self.samples[s].datamc == 0 : scale = 1.
#			print '%-12s: %f' % (s, self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_nloose_'+ratiovar).Integral()*scale)
			if i is 0 :
				h_ntight = self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_ntight_'+ratiovar).Clone()
				h_nloose = self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_nloose_'+ratiovar).Clone()
				h_ntight.Scale(scale)
				h_nloose.Scale(scale)
			else :
				h_ntight.Add(self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_ntight_'+ratiovar), scale)
				h_nloose.Add(self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_nloose_'+ratiovar), scale)
		return (h_ntight, h_nloose)


	def fit_fRatioPlots(self, h_data, h_wjets, h_zjets, h_qcd, individual_fit = True) :
		var = ROOT.RooRealVar('var', 'var', 0., 200.)

		h_ewk = h_wjets.Clone()
		h_ewk.Add(h_zjets)

		rdh_data  = ROOT.RooDataHist('data' ,'data' , ROOT.RooArgList(var), h_data , 1.)
		rdh_wjets = ROOT.RooDataHist('wjets','wjets', ROOT.RooArgList(var), h_wjets)
		rdh_zjets = ROOT.RooDataHist('zjets','zjets', ROOT.RooArgList(var), h_zjets)
		rdh_ewk   = ROOT.RooDataHist('ewk'  ,'ewk'  , ROOT.RooArgList(var), h_ewk  )
		rdh_qcd   = ROOT.RooDataHist('qcd'  ,'qcd'  , ROOT.RooArgList(var), h_qcd  )

		pdf_wjets = ROOT.RooHistPdf('wjets_pdf', 'wjets_pdf', ROOT.RooArgSet(var), rdh_wjets)
		pdf_zjets = ROOT.RooHistPdf('zjets_pdf', 'zjets_pdf', ROOT.RooArgSet(var), rdh_zjets)
		pdf_ewk   = ROOT.RooHistPdf('ewk_pdf'  , 'ewk_pdf'  , ROOT.RooArgSet(var), rdh_ewk  )
		pdf_qcd   = ROOT.RooHistPdf('qcd_pdf'  , 'qcd_pdf'  , ROOT.RooArgSet(var), rdh_qcd  )

		init_wjets = h_wjets.Integral()
		init_zjets = h_zjets.Integral()
		init_ewk   = h_ewk  .Integral()
		init_qcd   = h_qcd  .Integral()

		n_wjets = ROOT.RooRealVar('n_wjets', 'number of wjets', init_wjets, 0.70*init_wjets, 1.30*init_wjets)
		n_zjets = ROOT.RooRealVar('n_zjets', 'number of zjets', init_zjets, 0.70*init_zjets, 1.30*init_zjets)
		n_ewk   = ROOT.RooRealVar('n_ewk'  , 'number of ewk'  , init_ewk  , 0.50*init_ewk  , 2.00*init_ewk  )
		n_qcd   = ROOT.RooRealVar('n_qcd'  , 'number of qcd'  , init_qcd  , 0.50*init_qcd  , 2.00*init_qcd  )

		if individual_fit :
			model = ROOT.RooAddPdf('model', 'model', ROOT.RooArgList(pdf_wjets, pdf_zjets, pdf_qcd), ROOT.RooArgList(n_wjets, n_zjets, n_qcd))
		else :
			model = ROOT.RooAddPdf('model', 'model', ROOT.RooArgList(pdf_ewk, pdf_qcd), ROOT.RooArgList(n_ewk, n_qcd))
		model.fitTo(rdh_data, ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.Extended(), ROOT.RooFit.PrintLevel(-1))

		scale = {}
		if individual_fit :
			scale['wjets'] = n_wjets.getVal() / init_wjets
			scale['zjets'] = n_zjets.getVal() / init_zjets
		else :
			scale['wjets'] = n_ewk  .getVal() / init_ewk
			scale['zjets'] = scale['wjets']
		scale['qcd'  ] = n_qcd  .getVal() / init_qcd

		return scale

		mc = ROOT.TObjArray(3)
		mc.Add(h_wjets)
		mc.Add(h_zjets)
		mc.Add(h_qcd  )

		fit = ROOT.TFractionFitter(h_data, mc)
#		fit.Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
#		fit.SetRangeX(1,15);                    // use only the first 15 bins in the fit
		status = fit.Fit()
		print 'fit status:', status


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


	def plot_ratios(self, subdir = 'Ratios') :
		print '[status] producing fake and prompt ratio plots..'
#			(self.h2_MufRatio   , self.h_MufRatio_pt   , self.h_MufRatio_eta   , self.h_MufRatio_nv   , self.MufRatio   , self.MufRatioE   ) = self.calculateRatio(mu_samples, 'MM', 'SigSup', applyEwkSubtr , EWK_SF)
#			(self.h2_ElfRatio   , self.h_ElfRatio_pt   , self.h_ElfRatio_eta   , self.h_ElfRatio_nv   , self.ElfRatio   , self.ElfRatioE   ) = self.calculateRatio(el_samples, 'EE', 'SigSup', applyEwkSubtr , EWK_SF)
#			(self.h2_MupRatio   , self.h_MupRatio_pt   , self.h_MupRatio_eta   , self.h_MupRatio_nv   , self.MupRatio   , self.MupRatioE   ) = self.calculateRatio(mu_samples, 'MM', 'ZDecay')
#			(self.h2_ElpRatio   , self.h_ElpRatio_pt   , self.h_ElpRatio_eta   , self.h_ElpRatio_nv   , self.ElpRatio   , self.ElpRatioE   ) = self.calculateRatio(el_samples, 'EE', 'ZDecay')

		for TeX_switch in [True, False] :
			cms_label = 0
			pl = ttvStyle.ttvStyle(cms_label = cms_label, TeX_switch = TeX_switch)
			pl.lumi = self.lumi_HLTMu24Eta2p1
			self.save_ratio_plot(pl, '%s%s' % (self.path, subdir), self.h_MufRatio_pt , self.h_MufRatio_pt_MC , 'MufRatio_pt' , pl.get_varName('pTm'  ), pl.get_varName('fratio'))
			self.save_ratio_plot(pl, '%s%s' % (self.path, subdir), self.h_MufRatio_eta, self.h_MufRatio_eta_MC, 'MufRatio_eta', pl.get_varName('etam' ), pl.get_varName('fratio'))
			self.save_ratio_plot(pl, '%s%s' % (self.path, subdir), self.h_MufRatio_nv , self.h_MufRatio_nv_MC , 'MufRatio_nv' , pl.get_varName('NVrtx'), pl.get_varName('fratio'))
			pl.lumi = self.lumi
			self.save_ratio_plot(pl, '%s%s' % (self.path, subdir), self.h_MupRatio_pt , self.h_MupRatio_pt_MC , 'MupRatio_pt' , pl.get_varName('pTm'  ), pl.get_varName('pratio'))
			self.save_ratio_plot(pl, '%s%s' % (self.path, subdir), self.h_MupRatio_eta, self.h_MupRatio_eta_MC, 'MupRatio_eta', pl.get_varName('etam' ), pl.get_varName('pratio'))
			self.save_ratio_plot(pl, '%s%s' % (self.path, subdir), self.h_MupRatio_nv , self.h_MupRatio_nv_MC , 'MupRatio_nv' , pl.get_varName('NVrtx'), pl.get_varName('pratio'))
			pl.lumi = self.lumi_HLTEl17Jet30
			self.save_ratio_plot(pl, '%s%s' % (self.path, subdir), self.h_ElfRatio_pt , self.h_ElfRatio_pt_MC , 'ElfRatio_pt' , pl.get_varName('pTe'  ), pl.get_varName('fratio'))
			self.save_ratio_plot(pl, '%s%s' % (self.path, subdir), self.h_ElfRatio_eta, self.h_ElfRatio_eta_MC, 'ElfRatio_eta', pl.get_varName('etae' ), pl.get_varName('fratio'))
			self.save_ratio_plot(pl, '%s%s' % (self.path, subdir), self.h_ElfRatio_nv , self.h_ElfRatio_nv_MC , 'ElfRatio_nv' , pl.get_varName('NVrtx'), pl.get_varName('fratio'))
			pl.lumi = self.lumi
			self.save_ratio_plot(pl, '%s%s' % (self.path, subdir), self.h_ElpRatio_pt , self.h_ElpRatio_pt_MC , 'ElpRatio_pt' , pl.get_varName('pTe'  ), pl.get_varName('pratio'))
			self.save_ratio_plot(pl, '%s%s' % (self.path, subdir), self.h_ElpRatio_eta, self.h_ElpRatio_eta_MC, 'ElpRatio_eta', pl.get_varName('etae' ), pl.get_varName('pratio'))
			self.save_ratio_plot(pl, '%s%s' % (self.path, subdir), self.h_ElpRatio_nv , self.h_ElpRatio_nv_MC , 'ElpRatio_nv' , pl.get_varName('NVrtx'), pl.get_varName('pratio'))
		#pl.save_plot_2d(self.h2_MufRatio, 'p_{T} [GeV]', '#eta')


	def save_ratio_plot(self, pl, path, h_data, h_mc, name = '', x_title = '', y_title = '', x_range = None) :
		if name == '' : name = h_mc.GetName()
		h_data = h_data.Clone()
		h_mc   = h_mc  .Clone()
		h_data.UseCurrentStyle()
		h_mc  .UseCurrentStyle()
		if not x_range is None :
			h_data.SetAxisRange(x_range[0], x_range[1], 'X')
			h_mc  .SetAxisRange(x_range[0], x_range[1], 'X')
		h_mc.SetMarkerStyle(23)
		h_mc.SetMarkerColor(ROOT.kRed)
		h_mc.SetLineColor(ROOT.kRed)
		h_data.SetMinimum(0.)
		h_mc  .SetMinimum(0.)
		pl.get_maximum([h_data, h_mc], scale = 1.4, set_maximum = True)
		canvas = pl.get_canvas()
		canvas.cd()
		leg_entries = [[h_data, 'Data', 'lp'], [h_mc, 'Simulation', 'lp']]
		leg = pl.draw_legend(leg_entries)
		h_mc.Draw()
		h_mc.GetXaxis().SetTitle(x_title)
		h_mc.GetYaxis().SetTitle(y_title)
		h_data.Draw('same pe')
		pl.draw_cmsLine()
		leg.Draw()
		pl.save_canvas(canvas, path, name)


	def make_controlPlots(self, samples_data, ch_str) :
		variables = []
#		variables.append('NJets'      )
#		variables.append('HT'         )
#		variables.append('MaxJPt'     )
#		variables.append('NVertices'  )
#		variables.append('ClosJetPt'  )
#		variables.append('AwayJetPt'  )
#		variables.append('NBJets'     )
		variables.append('MET'        )
#		variables.append('MT'         )
#		variables.append('MET_noMTCut')
		variables.append('MT_MET30'   )
#		variables.append('LepPt'      )
		variables.append('LepEta'     )
#		variables.append('LepIso'     )
#		variables.append('ClosJetDR'  )
#		variables.append('AwayJetDR'  )

		ewk_sf = self.get_EWK_SF(samples_data, ch_str, do_fit = True)

		for var in variables :
			self.make_controlPlot(samples_data, ch_str, var, 1.)
			self.make_controlPlot(samples_data, ch_str, var, ewk_sf)


	def make_controlPlot(self, samples_data, chan_str, ratiovar, EWK_SF = 1.) :
		scale_str = ''
		if EWK_SF == 1. : scale_str = 'unscaled/'
		if type(EWK_SF) is dict : scale_str = 'fitted/'
		subdir = 'RatioControlPlots/%s%s' % (scale_str, chan_str)
		path = '%s%s/' % (self.path, subdir)
		helper.mkdir(path)
		if chan_str == 'el'   : lumi = self.lumi_HLTEl17Jet30
		if chan_str == 'mu17' : lumi = self.lumi_HLTMu17
		if chan_str == 'mu24' : lumi = self.lumi_HLTMu24Eta2p1

		processes = []
		processes.append('wjets')
		processes.append('zjets')
		processes.append('qcd'  )
		processes.append('data' )

		samples_wjets = sample.sample.get_samples('WJets' , self.samples)
		samples_zjets = sample.sample.get_samples('DYJets', self.samples)
		samples_qcd   = sample.sample.get_samples('QCD'   , self.samples)
#		samples_zjets = ['DYJets']

		histos = {}
		histos['tight'] = {}
		histos['loose'] = {}
		hstack = {}

		(histos['tight']['data' ], histos['loose']['data' ]) = self.get_fRatioPlots(samples_data , chan_str, ratiovar)
		(histos['tight']['wjets'], histos['loose']['wjets']) = self.get_fRatioPlots(samples_wjets, chan_str, ratiovar)
		(histos['tight']['zjets'], histos['loose']['zjets']) = self.get_fRatioPlots(samples_zjets, chan_str, ratiovar)
		(histos['tight']['qcd'  ], histos['loose']['qcd'  ]) = self.get_fRatioPlots(samples_qcd  , chan_str, ratiovar)

		for tl in histos :
			if type(EWK_SF) is dict :
				histos[tl]['wjets'].Scale(EWK_SF[tl]['wjets'])
				histos[tl]['zjets'].Scale(EWK_SF[tl]['zjets'])
				histos[tl]['qcd'  ].Scale(EWK_SF[tl]['qcd'  ])
			else :
				histos[tl]['wjets'].Scale(EWK_SF)
				histos[tl]['zjets'].Scale(EWK_SF)
			maximum = ttvStyle.ttvStyle.get_maximum(histos[tl].values())
			for TeX_switch in [True, False] :
				pl = ttvStyle.ttvStyle(lumi = lumi, cms_label = 0, TeX_switch = TeX_switch, short_names = False)
				canvas = pl.get_canvas()

				hstack = ROOT.THStack('hs_%s' % ratiovar, '%s' % ratiovar)
				leg_entries = []
				for process in processes :
					histos[tl][process].SetFillColor(pl.get_fillColor(process))
					if not (process == 'data' or process == 'obs') :
						hstack.Add(histos[tl][process])
				for process in reversed(processes) :
					if process == 'data' or process == 'obs' :
						leg_entries.append([histos[tl][process], pl.get_processName('obs'), 'lp'])
					else :
						leg_entries.append([histos[tl][process], pl.get_processName(process), 'f'])
				hstack.SetMaximum(maximum)
				leg = pl.draw_legend(leg_entries)
				hstack.Draw('hist')
				bin_width = hstack.GetXaxis().GetBinWidth(1)
				if ratiovar == 'NJets' or ratiovar == 'NBJets' or hstack.GetXaxis().IsVariableBinSize() :
					y_title = 'Events'
				else :
					y_title = pl.get_eventsPerGeVString(bin_width)
				hstack.GetXaxis().SetTitle(pl.get_varName(ratiovar))
				hstack.GetYaxis().SetTitle(y_title)
				histos[tl]['data'].Draw('psame')
				pl.draw_cmsLine()
				canvas.cd()
				leg.Draw()
				canvas.UseCurrentStyle()
				pl.save_canvas(canvas, path, 'N%s_%s' % (tl, ratiovar))
