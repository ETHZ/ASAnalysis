#!/usr/bin/python

import ROOT


def get_NJets(event) :
	njets = 0
	for i in range(event.NJets) :
		if is_GoodJet(event, i) : njets += 1
	return njets

def is_GoodJet(event, jet) :
	if event.JetPt[jet] < 30. : return False
	return True

def is_bTagged(event, jet) :
	if is_GoodJet(event, jet) and event.JetCSVBTag[jet] > 0.679 : return True
	return False


file = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/lbaeni/SSDLTrees/2013/Sep08/MC/TT-CT10-TuneZ2star-8TeV-powheg-tauola-Summer12-DR53X-PU-S10-START53-V7A-v1/output_0.root")
tree = file.Get("Analysis")

h_mjj = ROOT.TH1F("mjj", "2 jet mass", 50, 0., 500.)
h_Wcount  = ROOT.TH1F('Wcout', 'number djets with in the W mass', 10, 0., 10.)
h_bJcount = ROOT.TH1F('bJcount', 'number of b-jets', 10, 0., 10.)

counter = 0
#nevents = tree.GetEntries()

for event in tree :
	if counter > 10000 : break
	counter += 1
	W_counter = 0
	jetpair_counter = 0
	bJet_counter = 0
	dijetmass = 0.
	ind_jet1 = -1
	ind_jet2 = -1
#	if 100*counter/nevents%10 is 0. : print '[status]', 'done with %3.1f \% of events.' %(100.*counter/nevents)
	njets = get_NJets(event)
	if njets < 2 : continue
	for i in range(event.NJets) :
		if not is_GoodJet(event, i) : continue
		if is_bTagged(event, i) : continue
		for j in range(i+1, event.NJets) :
			if not is_GoodJet(event, j) : continue
			if is_bTagged(event, j) : continue
			jetpair_counter += 1
			jet1  = ROOT.TLorentzVector()
			jet2  = ROOT.TLorentzVector()
			dijet = ROOT.TLorentzVector()
			Wcand = ROOT.TLorentzVector()
			jet1.SetPtEtaPhiE(event.JetPt[i], event.JetEta[i], event.JetPhi[i], event.JetEnergy[i])
			jet2.SetPtEtaPhiE(event.JetPt[j], event.JetEta[j], event.JetPhi[j], event.JetEnergy[j])
			dijet = jet1 + jet2
			h_mjj.Fill(dijet.M())
			if dijet.M() < 5. :
				print '-------'
				print 'jet1:', event.JetPt[i], event.JetEta[i], event.JetPhi[i], event.JetEnergy[i]
				print 'jet2:', event.JetPt[j], event.JetEta[j], event.JetPhi[j], event.JetEnergy[j]
			if abs(dijet.M() - 80.) < abs(dijetmass - 80.) or dijetmass is 0. :
				dijetmass = dijet.M()
				ind_jet1 = i
				ind_jet2 = j
			if dijet.M() > 65. and dijet.M() < 95. :
				W_counter += 1
				Wcand = dijet
			if counter < 10 : print dijet.M()

	for i in range(event.NJets) :
		if not is_bTagged(event, i) : continue
		if i in (ind_jet1, ind_jet2) : continue # sanity check
		bJet_counter += 1
	h_bJcount.Fill(bJet_counter)

#	if dijetmass < 5. : print dijetmass, njets
#	h_mjj.Fill(dijetmass)
#	h_Wcount.Fill(W_counter)
	h_Wcount.Fill(jetpair_counter)

c1 = ROOT.TCanvas("canvas", "canvas", 0, 0, 1500, 1000)
c1.Divide(3, 2)
c1.cd(1)
h_mjj.Draw()
c1.cd(2)
h_Wcount.Scale(1./h_Wcount.Integral())
h_Wcount.Draw()
c1.cd(3)
h_bJcount.Draw()

raw_input('ok? ')



#		if counter < 10 : print event.Event, i, event.JetPt[i]




#float SSDLDumper::getM3(){
#	// Return M3
#	int njets = getNJets();
#	if(njets < 3) return -1.;
#
#	float triJetPtMax = 0.;
#	float m3 = 0.;
#	for(size_t i = 0; i < NJets; ++i){
#		if(!isGoodJet(i)) continue;
#		for(size_t j = i+1; j < NJets; ++j){
#			if(!isGoodJet(j)) continue;
#			for(size_t k = j+1; k < NJets; ++k){
#				if(!isGoodJet(k)) continue;
#				TLorentzVector jet1; jet1.SetPtEtaPhiE(JetPt[i], JetEta[i], JetPhi[i], JetEnergy[i]);
#				TLorentzVector jet2; jet2.SetPtEtaPhiE(JetPt[j], JetEta[j], JetPhi[j], JetEnergy[j]);
#				TLorentzVector jet3; jet3.SetPtEtaPhiE(JetPt[k], JetEta[k], JetPhi[k], JetEnergy[k]);
#				TLorentzVector triJet = jet1 + jet2 + jet3;
#				if( triJet.Pt() > triJetPtMax ){
#					triJetPtMax = triJet.Pt();
#					m3 = triJet.M();
#				}
#			}
#		}
#	}	
#	return m3;
#}
