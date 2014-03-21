#! /usr/bin/python
#import os, sys, commands, subprocess, math
import ROOT
import sys
import selection as sel
import sample

class plotter :
	'''read sigtree and produces plots'''

	def __init__(self, path) :
		print '[status] initialize plotter..'
		self.path = path
		self.ssdlfile = ROOT.TFile.Open(path + '/SSDLYields.root', 'READ')
		self.sigtree = self.ssdlfile.Get('SigEvents')
		print '[status] loaded SigEventsTree with %d events' % (self.sigtree.GetEntries())

		ROOT.gSystem.Load('./FakeRatios.so')

		# selections
		presel = sel.selection(name = 'presel', minNjets = 3, minNbjetsL = 1, minNbjetsM = 1)
		#print presel.get_selectionString()
		self.selections = {}
		self.selections[presel.name] = presel

		# samples
		self.samples = {}

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

		el_EWK_SF   = self.get_EWK_SF('el')
		mu17_EWK_SF = self.get_EWK_SF('mu17')
		mu24_EWK_SF = self.get_EWK_SF('mu24')

		#(h2_tight, h2_loose) = self.get_TightLoose(self.get_samples('DoubleEle'), 'EE', 'SigSup')
		(h2_tight, h2_loose) = self.get_TightLoose(['DYJets'], 'MM', 'ZDecay')
		h2_tight.Draw('box')
		raw_input('ok? ')

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
		return samplelist


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
				h_ntight = self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_ntight_'+ratiovar)
				h_nloose = self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_nloose_'+ratiovar)
				h_ntight.Scale(scale)
				h_nloose.Scale(scale)
			else :
				h_ntight.Add(self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_ntight_'+ratiovar), scale)
				h_nloose.Add(self.ssdlfile.Get(s+'/FRatioPlots/'+s+'_'+chan_str+'_nloose_'+ratiovar), scale)
		return (h_ntight, h_nloose)


	def get_ChMisID_SF(self) :
		foo = 0


	def fill_ratios(self) :
		print '[status] filling fake and prompt ratio histograms..'
		## void SSDLPlotter::fillRatios(vector<int> musamples, vector<int> elsamples,
		## 			     int datamc, bool applyEwkSubtr, bool printOutput){
		## 	if(datamc == 0){
		## 	  cout << "Filling ratios for Data... " << endl; 
		## NOT NEEDED 	  fH1D_MufRatio = fillRatioPt(Muon, musamples, SigSup, applyEwkSubtr, printOutput);
		## NOT NEEDED 	  fH1D_MupRatio = fillRatioPt(Muon, musamples, ZDecay, applyEwkSubtr, printOutput);
		## NOT NEEDED 	  fH1D_ElfRatio = fillRatioPt(Elec, elsamples, SigSup, applyEwkSubtr, printOutput);
		## NOT NEEDED 	  fH1D_ElpRatio = fillRatioPt(Elec, elsamples, ZDecay, applyEwkSubtr, printOutput);
		self.fH2D_MufRatio = fillRatio('Muon', musamples, 'SigSup', applyEwkSubtr, printOutput)
		self.fH2D_MupRatio = fillRatio('Muon', musamples, 'ZDecay', applyEwkSubtr, printOutput)
		self.fH2D_ElfRatio = fillRatio('Elec', self.get_samples('DoubleEle'), 'SigSup', applyEwkSubtr, printOutput)
		self.fH2D_ElpRatio = fillRatio('Elec', self.get_samples('DoubleEle'), 'ZDecay', applyEwkSubtr, printOutput)
		## 	}
		## NOT NEEDED 	if(datamc == 1){
		## NOT NEEDED 	  cout << "Filling ratios for MC... " << endl;
		## NOT NEEDED 	  fH1D_MufRatio_MC = fillRatioPt(Muon, musamples, SigSup, applyEwkSubtr, printOutput);
		## NOT NEEDED 	  fH1D_MupRatio_MC = fillRatioPt(Muon, musamples, ZDecay, applyEwkSubtr, printOutput);
		## NOT NEEDED 	  fH1D_ElfRatio_MC = fillRatioPt(Elec, elsamples, SigSup, applyEwkSubtr, printOutput);
		## NOT NEEDED 	  fH1D_ElpRatio_MC = fillRatioPt(Elec, elsamples, ZDecay, applyEwkSubtr, printOutput);
		## NOT NEEDED 	  fH2D_MufRatio_MC = fillRatio(  Muon, musamples, SigSup, applyEwkSubtr, printOutput);
		## NOT NEEDED 	  fH2D_MupRatio_MC = fillRatio(  Muon, musamples, ZDecay, applyEwkSubtr, printOutput);
		## NOT NEEDED 	  fH2D_ElfRatio_MC = fillRatio(  Elec, elsamples, SigSup, applyEwkSubtr, printOutput);
		## NOT NEEDED 	  fH2D_ElpRatio_MC = fillRatio(  Elec, elsamples, ZDecay, applyEwkSubtr, printOutput);
		## NOT NEEDED 	}
		## }


	def fill_ratio(self, chan, samples, region, applyEwkSubtr) :
		foo = 0
		# sets up the histograms and calls calculateRatio

		## TH2D* SSDLPlotter::fillRatio(gChannel chan, vector<int> samples, gFPSwitch fp, bool applyEwkSubtr, bool output){
		## 	gStyle->SetOptStat(0);
		## 	TString shortname[2] = {"Mu", "El"};
		## 	TString longname[2] = {"Muons", "Electrons"};
		## 	int muelswitch = 0;
		## 	if(chan == Elec) muelswitch = 1;
		## 
		## 	TH2D *h_2d;
		## 	TH1D *h_pt, *h_eta;
		## 	if(fp == SigSup){
		## 		h_2d  = new TH2D(Form("%sRatio",   shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt vs Eta", longname[muelswitch].Data()), getNFPtBins(chan), getFPtBins(chan), getNEtaBins(chan), getEtaBins(chan));
		## 		h_pt  = new TH1D(Form("%sRatioPt" ,shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt"       , longname[muelswitch].Data()), getNFPtBins(chan), getFPtBins(chan));
		## 		h_eta = new TH1D(Form("%sRatioEta",shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Eta"      , longname[muelswitch].Data()), getNEtaBins(chan), getEtaBins(chan));
		## 	}
		## 	if(fp == ZDecay){
		## 		h_2d  = new TH2D(Form("%sRatio",   shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt vs Eta", longname[muelswitch].Data()), getNPPtBins(chan), getPPtBins(chan), getNEtaBins(chan), getEtaBins(chan));
		## 		h_pt  = new TH1D(Form("%sRatioPt" ,shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Pt"       , longname[muelswitch].Data()), getNPPtBins(chan), getPPtBins(chan));
		## 		h_eta = new TH1D(Form("%sRatioEta",shortname[muelswitch].Data()), Form("Ratio of tight to loose %s vs Eta"      , longname[muelswitch].Data()), getNEtaBins(chan), getEtaBins(chan));
		## 	}
		## 
		## 	h_2d->SetXTitle("p_{#perp} (GeV)");
		## 	h_2d->SetYTitle("#eta");
		## 	h_2d->SetZTitle("# Tight / # Loose");
		## 
		## 	calculateRatio(samples, chan, fp, h_2d, h_pt, h_eta, applyEwkSubtr, output);
		## 	delete h_pt, h_eta;
		## 	return h_2d;
		## }


	def calculateRatio(self, applyEwkSubtr, EWK_SF) :
		foo = 0
		# - sets up ratio histos
		# - calls getPassedTotal / get_TightLoose to get ntight and nloose histos

		## void SSDLPlotter::calculateRatio(vector<int> samples, gChannel chan, gFPSwitch fp, 
		## 				 TH2D*& h_2d, TH1D*& h_pt, TH1D*& h_eta, TH1D*& h_nv, 
		## 				 bool applyEwkSubtr, bool output, bool ttbarMatched){
		## /*
		## TODO Fix treatment of statistical errors and luminosity scaling here!
		## */
		## 	gStyle->SetOptStat(0);
		## 
		## 	h_2d->Sumw2();
		## 	h_pt->Sumw2();
		## 	h_eta->Sumw2();
		## 
		## 	vector<int> wjets_samples, zjets_samples;
		## 	wjets_samples.push_back(WJets);
		## 	zjets_samples.push_back(DYJets);
		## 
		## 	TH2D *H_ntight = new TH2D("NTight", "NTight Muons", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
		## 	TH2D *H_nloose = new TH2D("NLoose", "NLoose Muons", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
		## 	H_ntight->Sumw2();
		## 	H_nloose->Sumw2();
		## 	TH2D *H_ntight_wjets = new TH2D("NTight_WJets", "NTight Muons WJets", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
		## 	TH2D *H_nloose_wjets = new TH2D("NLoose_WJets", "NLoose Muons WJets", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
		## 	H_ntight_wjets->Sumw2();
		## 	H_nloose_wjets->Sumw2();
		## 	TH2D *H_ntight_zjets = new TH2D("NTight_ZJets", "NTight Muons ZJets", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
		## 	TH2D *H_nloose_zjets = new TH2D("NLoose_ZJets", "NLoose Muons ZJets", h_2d->GetNbinsX(), h_2d->GetXaxis()->GetXbins()->GetArray(), h_2d->GetNbinsY(),  h_2d->GetYaxis()->GetXbins()->GetArray());
		## 	H_ntight_zjets->Sumw2();
		## 	H_nloose_zjets->Sumw2();
		## 
		## 	TH1D *H_ntight_nv = new TH1D("NTight_NV", "NTight Muons NV", 18, 0., 36.); H_ntight_nv->Sumw2();
		## 	TH1D *H_nloose_nv = new TH1D("NLoose_NV", "NLoose Muons NV", 18, 0., 36.); H_nloose_nv->Sumw2();
		## 	TH1D *H_ntight_nv_wjets = new TH1D("NTight_NV_WJets", "NTight Muons NV WJets", 18, 0., 36.); H_ntight_nv_wjets->Sumw2();
		## 	TH1D *H_nloose_nv_wjets = new TH1D("NLoose_NV_WJets", "NLoose Muons NV WJets", 18, 0., 36.); H_nloose_nv_wjets->Sumw2();
		## 	TH1D *H_ntight_nv_zjets = new TH1D("NTight_NV_ZJets", "NTight Muons NV ZJets", 18, 0., 36.); H_ntight_nv_zjets->Sumw2();
		## 	TH1D *H_nloose_nv_zjets = new TH1D("NLoose_NV_ZJets", "NLoose Muons NV ZJets", 18, 0., 36.); H_nloose_nv_zjets->Sumw2();
		## 	if (gRatiosFromTTbar || ttbarMatched) getPassedTotalTTbar(samples, chan, fp, H_ntight, H_nloose, output);
		## 	// marc jan 27 else getPassedTotal(samples, chan, fp, H_ntight, H_nloose, output);
		## 	else {
		## 		getPassedTotal(samples, chan, fp, H_ntight, H_nloose, H_ntight_nv, H_nloose_nv, output);
		## 		getPassedTotal(wjets_samples, chan, fp, H_ntight_wjets, H_nloose_wjets, H_ntight_nv_wjets, H_nloose_nv_wjets, output);
		## 		getPassedTotal(zjets_samples, chan, fp, H_ntight_zjets, H_nloose_zjets, H_ntight_nv_zjets, H_nloose_nv_zjets, output);
		## 	}
		## 	//if (fp == SigSup && gEWKCorrection) {
		## 	if (fp == SigSup && applyEwkSubtr) {
		## 		float lumi(1.);
		## 		//if (chan == Muon) lumi = fLumiNormHLTMu17      * fEWKMuSF; //Can we remove this line ?? BM
		## 		if (chan == Elec) lumi = fLumiNormHLTEl17Jet30 * fEWKElSF;
		## 		// mc samples are scaled to fLumiNorm. now scale to prescaled triggers.
		## 		float scale_vjets = lumi / fLumiNorm;
		## 		if (chan == Elec) {
		## 			H_ntight   ->Add(H_ntight_wjets   , (-1.) * scale_vjets);
		## 			H_ntight   ->Add(H_ntight_zjets   , (-1.) * scale_vjets);
		## 			H_nloose   ->Add(H_nloose_wjets   , (-1.) * scale_vjets);
		## 			H_nloose   ->Add(H_nloose_zjets   , (-1.) * scale_vjets);
		## 			H_ntight_nv->Add(H_ntight_nv_wjets, (-1.) * scale_vjets);
		## 			H_ntight_nv->Add(H_ntight_nv_zjets, (-1.) * scale_vjets);
		## 			H_nloose_nv->Add(H_nloose_nv_wjets, (-1.) * scale_vjets);
		## 			H_nloose_nv->Add(H_nloose_nv_zjets, (-1.) * scale_vjets);
		## 		}
		## 		if (chan == Muon) {
		## 			for (int ptbin = 1; ptbin < gNMuFPtBins+1; ptbin++) {
		## 				for (int etabin = 1; etabin < gNMuEtabins+1; etabin++) {
		## 					int bin = H_nloose->GetBin(ptbin, etabin);
		## 					//if (ptbin > 1 && etabin < 3) lumi = fLumiNormHLTMu24Eta2p1 * fEWKMu24SF;
		## 					if(H_nloose->GetXaxis()->GetBinLowEdge(ptbin)>=25.0 && 
		## 					   H_nloose->GetYaxis()->GetBinUpEdge(etabin)<=2.1 ) lumi = fLumiNormHLTMu24Eta2p1 * fEWKMu24SF;
		## 					else                         lumi = fLumiNormHLTMu17       * fEWKMu17SF;
		## 					float scale = lumi / fLumiNorm;
		## 					H_ntight   ->AddBinContent(bin, (-1.) * scale * H_ntight_wjets->GetBinContent(bin));
		## 					H_ntight   ->AddBinContent(bin, (-1.) * scale * H_ntight_zjets->GetBinContent(bin));
		## 					H_nloose   ->AddBinContent(bin, (-1.) * scale * H_nloose_wjets->GetBinContent(bin));
		## 					H_nloose   ->AddBinContent(bin, (-1.) * scale * H_nloose_zjets->GetBinContent(bin));
		## 				}
		## 			}
		## 		}
		## 	}
		## 	h_2d->Divide(H_ntight,    H_nloose,    1., 1., "B");
		## 	h_nv->Divide(H_ntight_nv, H_nloose_nv, 1., 1., "B");
		## 
		## 	TH1D *hmuloosept  = H_nloose->ProjectionX();
		## 	TH1D *hmulooseeta = H_nloose->ProjectionY();
		## 	TH1D *hmutightpt  = H_ntight->ProjectionX();
		## 	TH1D *hmutighteta = H_ntight->ProjectionY();
		## 
		## 	h_pt ->Divide(hmutightpt,  hmuloosept,  1., 1., "B"); // binomial
		## 	h_eta->Divide(hmutighteta, hmulooseeta, 1., 1., "B"); // weights are ignored
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
		## 	delete H_ntight, H_nloose, H_ntight_nv, H_nloose_nv, hmuloosept, hmulooseeta, hmutightpt, hmutighteta;
		## }


	def get_TightLoose(self, samples, chan_str, fp) :
		## chan_str: 'MM', 'EE'

		for i, s in enumerate(samples) :
			scale = 1.
			if self.samples[s].datamc != 0 : scale = self.lumi / self.samples[s].getLumi()
			if fp is 'SigSup' :
				h2_ntight   = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_fNTight')
				h2_nloose   = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_fNLoose')
				h_ntight_nv = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_fNTight_nv')
				h_nloose_nv = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_fNLoose_nv')
			elif fp is 'ZDecay' :
				h2_ntight   = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_pNTight')
				h2_nloose   = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_pNLoose')
				h_ntight_nv = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_pNTight_nv')
				h_nloose_nv = self.ssdlfile.Get(s+'/TLRatios/'+s+'_'+chan_str+'_pNLoose_nv')
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

		return (h2_passed, h2_total)
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


	def make_IntPredictions(self, selections) :
		'''oberservation and prediction for different selections'''
		foo = 0


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

	if ('-d' in args) and (args[args.index('-d')+1] is not '') :
		path = str(args[args.index('-d')+1])
		print path

	if ('-c' in args) and (args[args.index('-c')+1] is not '') :
		cardfile = str(args[args.index('-c')+1])
		print cardfile

	plotter = plotter(path)
	plotter.do_analysis(cardfile)
