#! /usr/bin/python
import ROOT
import os, sys
import helper


class ttvplot :

	def __init__(self, path, chan, lumi = 19500., cms_label = 0, asymmErr = True, TeX_switch = False) :
		self.path = path
		if not os.path.exists(self.path) :
			os.makedirs(self.path)
		self.chan = chan
		self.lumi = lumi
		self.cms_label = cms_label
		self.asymmErr = asymmErr

		# colors
		self.colors = {}
		self.colors['obs'  ] = 1
		self.colors['fake' ] = 46
		self.colors['chmid'] = 49
		self.colors['rare' ] = 38
		self.colors['wz'   ] = 39
		self.colors['ttz'  ] = 42
		self.colors['ttw'  ] = 44
		self.colors['btag' ] = 31
		self.colors['zz'   ] = 30

		# process names
		self.process_names = {}
		self.process_names['obs'  ] = 'Observed'
		self.process_names['fake' ] = 'Non-prompt Lepton'
		self.process_names['chmid'] = 'Charge MisID'
		self.process_names['rare' ] = 'Irreducible'
		self.process_names['wz'   ] = 'WZ'
		self.process_names['ttz'  ] = 't#bar{t}Z'
		self.process_names['ttw'  ] = 't#bar{t}W'
		self.process_names['bgtot'] = 'Backgrounds'
		self.process_names['btag' ] = 'Non-top'
		self.process_names['zz'   ] = 'ZZ'
		if TeX_switch is True :
			self.process_names['ttz'  ] = '\\ttbar{}Z'
			self.process_names['ttw'  ] = '\\ttbar{}W'

		# variable names
		self.var_names = {}
		self.var_names['HT'    ] = 'H_{T} [GeV]'
		self.var_names['MET'   ] = 'Particle Flow E_{T}^{miss} [GeV]'
		self.var_names['NJ'    ] = 'Jet Multiplicity'
		self.var_names['NbJmed'] = 'b-Jet Multiplicity (CSVM)'
		self.var_names['pT1'   ] = 'Leading Lepton p_{T} [GeV]'
		self.var_names['pT2'   ] = 'Sub-Leading Lepton p_{T} [GeV]'
		self.var_names['Int'   ] = ''
		self.var_names['Mll'   ] = 'm_{ll} [GeV]'
		self.var_names['NVrtx' ] = 'N_{Vertices}'
		self.var_names['minMT' ] = 'M_{T} [GeV]'
		self.var_names['M3'    ] = 'm_{jjj} [GeV]'


	def get_fillColor(self, process) :
		if process in self.colors.keys() :
			return self.colors[process]
		return 0


	def get_varName(self, var) :
		if var in self.var_names.keys() :
			return self.var_names[var]
		return ''


	def get_processName(self, process) :
		if process in self.process_names.keys() :
			return self.process_names[process]
		return '?'


	def save_plot(self, histos, var, selname = '', prefix = '', charge_str = '') :
		'''save plot with observation and predictions'''

		if selname != '' : selname += '_'
		if prefix  != '' : prefix = '_' + prefix

		# data with asymmetric errors
		gr_obs = helper.getGraphPoissonErrors(histos['obs'])

		# adding predictions to stack
		hs_pred = ROOT.THStack('hs_pred', 'hs_pred')
		if self.chan == '2L' :
			hs_pred.Add(histos['fake' ])
			hs_pred.Add(histos['chmid'])
			hs_pred.Add(histos['rare' ])
			hs_pred.Add(histos['wz'   ])
			hs_pred.Add(histos['ttz'  ])
			hs_pred.Add(histos['ttw'  ])
		elif self.chan == '3L' :
			hs_pred.Add(histos['fake' ])
			hs_pred.Add(histos['btag' ])
			hs_pred.Add(histos['rare' ])
			hs_pred.Add(histos['ttz'  ])
		elif self.chan == '4L' :
			hs_pred.Add(histos['fake' ])
			hs_pred.Add(histos['zz'   ])
			hs_pred.Add(histos['rare' ])
			hs_pred.Add(histos['ttz'  ])

		# set minimum and maximum
		scale = 1.8
		if var == 'Int' : scale = 2.4
		maximum = scale * max(histos['obs'].GetMaximum(), histos['pred'].GetMaximum())
		gr_obs .SetMaximum(maximum)
		hs_pred.SetMaximum(maximum)
		gr_obs .SetMinimum(0.)
		hs_pred.SetMinimum(0.)

		for process, histo in histos.iteritems() :
			histo.SetMaximum(maximum)
			histo.SetMinimum(0.)
			histo.SetLineColor(1)
			histo.SetLineWidth(1)
			histo.SetFillColor(self.get_fillColor(process))

		# data histogram settings
		histos['obs'  ].SetMarkerStyle(20)
		histos['obs'  ].SetMarkerSize(1.1)
		histos['obs'  ].SetLineWidth(2)
		gr_obs         .SetMarkerStyle(20)
		gr_obs         .SetMarkerSize(1.1)
		gr_obs         .SetLineWidth(2)
		
		# special settings for total BG line and uncertainty
		histos['bgtot'].SetLineWidth(3)
		histos['bgtot'].SetLineColor(1)
#		histos['bgtot'].SetFillColor(12)
		histos['bgtot'].SetFillStyle(0)

		histos['pred' ].SetLineWidth(0)
		histos['pred' ].SetFillColor(12)
		histos['pred' ].SetFillStyle(3005)

		# legend
		leg = ROOT.TLegend()
		leg.AddEntry(histos['obs'  ], self.process_names['obs'  ], 'lp')
		leg.AddEntry(histos['fake' ], self.process_names['fake' ], 'f')
		if self.chan == '2L' : leg.AddEntry(histos['chmid'], self.process_names['chmid'], 'f')
		if self.chan == '3L' : leg.AddEntry(histos['btag' ], self.process_names['btag' ], 'f')
		if self.chan == '4L' : leg.AddEntry(histos['zz'   ], self.process_names['zz'   ], 'f')
		leg.AddEntry(histos['rare' ], self.process_names['rare' ], 'f')
		if self.chan == '2L' : leg.AddEntry(histos['wz'   ], self.process_names['wz'   ], 'f')
		leg.AddEntry(histos['ttz'  ], self.process_names['ttz'  ], 'f')
		if self.chan == '2L' : leg.AddEntry(histos['ttw'  ], self.process_names['ttw'  ], 'f')
		leg.AddEntry(histos['bgtot'], self.process_names['bgtot'], 'l')
		leg.AddEntry(histos['pred' ], 'BG uncertainty', 'fl')
		# set position
		width = 0.17
		x = 0.65
		y = 0.93
		leg.SetX1NDC(x)
		leg.SetX2NDC(x+width)
		leg.SetY1NDC(y-leg.GetNRows()*0.25*width)
		leg.SetY2NDC(y)
		leg.SetFillStyle(0)
		leg.SetTextFont(42)
		leg.SetTextSize(0.03)
		leg.SetBorderSize(0)
		leg.SetTextAlign(12)

		canvas = ROOT.TCanvas('C_ObsPred_'+selname+var, 'Observed vs Predicted', 0, 0, 600, 600)
		canvas.SetLeftMargin(0.12)
		canvas.SetRightMargin(0.04)
		canvas.SetTopMargin(0.04)
		canvas.SetBottomMargin(0.12)
		canvas.cd()

#		ROOT.gStyle.SetOptStat(0)
		ROOT.gStyle.SetOptTitle(0)
		ROOT.gStyle.SetEndErrorSize(0)  # set the size of the small line at the end of the error bars
		ROOT.gPad.SetTicks(1,1)

		# draw and set axis titles
		hs_pred.Draw()
		hs_pred.GetXaxis().SetTitle(self.get_varName(var))
		hs_pred.GetYaxis().SetTitle('Events')
		hs_pred.GetXaxis().SetTitleOffset(1.25)
		hs_pred.GetYaxis().SetTitleOffset(1.25)
		hs_pred.GetXaxis().SetTitleSize(0.046)
		hs_pred.GetYaxis().SetTitleSize(0.046)
		hs_pred.GetXaxis().SetLabelSize(0.04)
		hs_pred.GetYaxis().SetLabelSize(0.04)
		#hs_pred.GetXaxis().SetNdivisions(206)
		if var == 'Int' :
			for bin in range(1, hs_pred.GetXaxis().GetNbins()+1) :
				if bin == 1 : binlabel = '#mu'+charge_str+'#mu'+charge_str
				if bin == 2 : binlabel = 'e'+charge_str+'#mu'+charge_str
				if bin == 3 : binlabel = 'e'+charge_str+'e'+charge_str
				hs_pred.GetXaxis().SetBinLabel(bin, binlabel)
				hs_pred.GetXaxis().SetLabelSize(0.062)
		elif 'NJ' in var or 'NbJ' in var :
			for bin in range(1, hs_pred.GetXaxis().GetNbins()+1) :
				binlabel = '%d' % (hs_pred.GetXaxis().GetBinLowEdge(bin))
				hs_pred.GetXaxis().SetBinLabel(bin, binlabel)
				hs_pred.GetXaxis().SetLabelSize(0.062)
		elif var == 'CFChan' :
			for bin in range(1, hs_pred.GetXaxis().GetNbins()+1) :
				if bin < 4 : charge_str = '^{+}'
				else       : charge_str = '^{-}'
				if bin == 1 or bin == 4 : binlabel = '#mu'+charge_str+'#mu'+charge_str
				if bin == 2 or bin == 5 : binlabel = 'e'+charge_str+'#mu'+charge_str
				if bin == 3 or bin == 6 : binlabel = 'e'+charge_str+'e'+charge_str
				hs_pred.GetXaxis().SetBinLabel(bin, binlabel)
				hs_pred.GetXaxis().SetLabelSize(0.062)
		elif var == 'Charge' :
			for bin in range(1, hs_pred.GetXaxis().GetNbins()+1) :
				if bin == 1 : binlabel = '--'
				if bin == 2 : binlabel = '++'
				hs_pred.GetXaxis().SetBinLabel(bin, binlabel)
				hs_pred.GetXaxis().SetLabelSize(0.062)
		elif not hs_pred.GetXaxis().IsVariableBinSize() and var != 'NVrtx' :
			bin_width = hs_pred.GetXaxis().GetBinWidth(1)
			hs_pred.GetYaxis().SetTitle('Events / %.0f GeV' % bin_width)
		##		if (diffVarName[var] == "NJ" || diffVarName[var] == "NbJmed"){
		##			for(size_t i = 1; i <= nbins[var]; ++i) hs_pred[var]->GetXaxis()->SetBinLabel(i, Form("%d", (int)bins[var][i-1]));
		##			hs_pred[var]->GetXaxis()->SetLabelSize(0.07);
		##			hs_pred[var]->GetXaxis()->SetTitleSize(0.045);
		##			hs_pred[var]->GetXaxis()->SetTitleOffset(1.07);
		##		}
		##		if (diffVarName[var] == "Int") {
		##			for(size_t i = 1; i <= nbins[var]; ++i) {
		##				TString binlabel = "?";
		##				if (i == 1) binlabel = "e"  + chargeSignString + "e"  + chargeSignString;
		##				if (i == 2) binlabel = "#mu"+ chargeSignString + "#mu"+ chargeSignString;
		##				if (i == 3) binlabel = "e"  + chargeSignString + "#mu"+ chargeSignString;
		##				hs_pred[var]->GetXaxis()->SetBinLabel(i, binlabel);
		##			}
		##		}
		leg.Draw()
		histos['pred'].Draw('0 E2 same')
		histos['bgtot'].Draw('hist same')
		if self.asymmErr :
			gr_obs         .Draw('PE same')
		else :
			histos['obs'  ].Draw('PE X0 same')
#		self.drawTopLine(0.56, 0.8)
		self.draw_cmsLine()
#		raw_input('ok? ')

		canvas.Print('%sObsPred%s_%s.pdf' % (self.path, prefix, var))
		canvas.Print('%sObsPred%s_%s.png' % (self.path, prefix, var))
#		raw_input('ok? ')


	def drawTopLine(self, rightedge = 0.60, scale = 1., leftedge = 0.13) :
		latex = ROOT.TLatex()
		latex.SetNDC()
		latex.SetTextFont(62)
		latex.SetTextSize(scale * 0.05)
		latex.DrawLatex(leftedge, 0.92, 'CMS Preliminary')
		latex.SetTextFont(42)
		latex.SetTextSize(scale * 0.04)
		latex.DrawLatex(rightedge, 0.92, '#sqrt{s} = 8 TeV, L_{int} = %4.1f fb^{-1}' % (self.lumi/1000.))


	def draw_cmsLine(self) :
		latex = ROOT.TLatex()
		latex.SetNDC()
		latex.SetTextFont(62)
		latex.SetTextSize(0.04)
		latex.SetTextAlign(13)
		if self.cms_label == 0 : cms_str = 'CMS'
		if self.cms_label == 1 : cms_str = 'CMS Simulation'
		if self.cms_label == 2 : cms_str = 'CMS Preliminary'
		if self.cms_label == 3 : cms_str = 'CMS Simulation Preliminary'
		latex.DrawLatex(0.15, 0.93, cms_str)
		latex.SetTextFont(42)
		latex.SetTextSize(0.03)
		latex.DrawLatex(0.15, 0.88, '#sqrt{s} = 8 TeV, L_{int} = %4.1f fb^{-1}' % (self.lumi/1000.))


	def read_histos(self, path, var) :
		file = ROOT.TFile.Open(path, 'READ')
		histos = {}

		if self.chan == '2L' :
			histos['obs'  ] = file.Get('h_obs'  +'_'+var).Clone()
			histos['fake' ] = file.Get('h_fake' +'_'+var).Clone()
			histos['chmid'] = file.Get('h_chmid'+'_'+var).Clone()
			histos['rare' ] = file.Get('h_rare' +'_'+var).Clone()
			histos['wz'   ] = file.Get('h_wz'   +'_'+var).Clone()
			histos['ttz'  ] = file.Get('h_ttz'  +'_'+var).Clone()
			histos['ttw'  ] = file.Get('h_ttw'  +'_'+var).Clone()
			histos['bgtot'] = file.Get('h_bgtot'+'_'+var).Clone()
			histos['pred' ] = file.Get('h_pred' +'_'+var).Clone()

		if self.chan == '3L' :
			histos['obs'  ] = file.Get('h_obs'  +'_'+var).Clone()
			histos['fake' ] = file.Get('h_fake' +'_'+var).Clone()
			histos['btag' ] = file.Get('h_btag' +'_'+var).Clone()
			histos['rare' ] = file.Get('h_rare' +'_'+var).Clone()
			histos['ttz'  ] = file.Get('h_ttz'  +'_'+var).Clone()
			histos['bgtot'] = file.Get('h_bgtot'+'_'+var).Clone()
			histos['pred' ] = file.Get('h_pred' +'_'+var).Clone()

		if self.chan == '4L' :
			histos['obs'  ] = file.Get('h_obs'  +'_'+var).Clone()
			histos['fake' ] = file.Get('h_fake' +'_'+var).Clone()
			histos['zz'   ] = file.Get('h_zz'   +'_'+var).Clone()
			histos['rare' ] = file.Get('h_rare' +'_'+var).Clone()
			histos['ttz'  ] = file.Get('h_ttz'  +'_'+var).Clone()
			histos['bgtot'] = file.Get('h_bgtot'+'_'+var).Clone()
			histos['pred' ] = file.Get('h_pred' +'_'+var).Clone()

		return histos


if __name__ == '__main__' :
	args = sys.argv

	if ('--help' in args) or ('-h' in args) or ('-i' not in args) or ('-o' not in args) or ('--2L' not in args and '--3L' not in args and '--4L' not in args) :
		print 'usage: python -i <INPUT ROOTFILE> -o <OUTPUTDIR> --2L/--3L/--4L'
		sys.exit(1)

	if ('-i' in args) and (args[args.index('-i')+1] != '') :
		inputfile = args[args.index('-i')+1]

	if ('-o' in args) and (args[args.index('-o')+1] != '') :
		outputpath = args[args.index('-o')+1]

	if   '--2L' in args : chan = '2L'
	elif '--3L' in args : chan = '3L'
	elif '--4L' in args : chan = '4L'

	pl = ttvplot(outputpath, chan)
	histos = pl.read_histos(inputfile, 'Mll')
	pl.save_plot(histos, outputpath, 'Mll')
	# TODO: loop over all variables in the root file
