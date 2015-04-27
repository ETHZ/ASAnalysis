#! /usr/bin/python
import ROOT
import os, sys
import helper
import ttvStyle


class ttvplot(ttvStyle.ttvStyle) :

	def __init__(self, path, chan, lumi = -1., cms_label = 0, asymmErr = True, TeX_switch = False, short_names = False) :
		ttvStyle.ttvStyle.__init__(self, lumi = lumi, cms_label = cms_label, TeX_switch = TeX_switch, short_names = short_names)
		self.path = path
		if not self.path.endswith('/') : self.path += '/'
		if not os.path.exists(self.path) :
			os.makedirs(self.path)
		self.chan = chan
		self.asymmErr = asymmErr


	def save_plot(self, histos, var, prefix = '', suffix = '', charge_str = '') :
		'''save plot with observation and predictions'''

		# clone all histograms to avoid that histograms get scaled multiple times
		histos_tmp = {}
		for histo in histos :
			histos_tmp[histo] = histos[histo].Clone()
		histos = histos_tmp

		if prefix  != '' and not prefix.startswith('_') : prefix = '_' + prefix
		if suffix  != '' and not suffix.startswith('_') : suffix = '_' + suffix

		# data with asymmetric errors
#		gr_obs = helper.getGraphPoissonErrors_new(histos['obs'])

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

		# get canvas
		canvas = self.get_canvas('C_ObsPred')
		canvas.cd()
		self.ttvStyle.cd()

		# set styles
		for process, histo in histos.iteritems() :
			histo.SetFillColor(self.get_fillColor(process))

		# legend
		leg_entries = []
		leg_entries.append([histos['obs'  ], self.get_processName('obs'  ), 'lp'])
		if self.chan == '3L' :
			leg_entries.append([histos['ttz'  ], self.get_processName('ttz'  ), 'f'])
			leg_entries.append([histos['ttw'  ], self.get_processName('ttw'  ), 'f'])
			leg_entries.append([histos['rare' ], self.get_processName('rare' ), 'f'])
			leg_entries.append([histos['btag' ], self.get_processName('btag' ), 'f'])
		if self.chan == '4L' :
			leg_entries.append([histos['ttz'  ], self.get_processName('ttz'  ), 'f'])
			leg_entries.append([histos['ttw'  ], self.get_processName('ttw'  ), 'f'])
			leg_entries.append([histos['rare' ], self.get_processName('rare' ), 'f'])
			leg_entries.append([histos['zz'   ], self.get_processName('zz'   ), 'f'])
		if self.chan == '2L' :
			leg_entries.append([histos['ttw'  ], self.get_processName('ttw'  ), 'f'])
			leg_entries.append([histos['ttz'  ], self.get_processName('ttz'  ), 'f'])
			leg_entries.append([histos['wz'   ], self.get_processName('wz'   ), 'f'])
			leg_entries.append([histos['rare' ], self.get_processName('rare' ), 'f'])
			leg_entries.append([histos['chmid'], self.get_processName('chmid'), 'f'])
		leg_entries.append([histos['fake' ], self.get_processName('fake' ), 'f'])
		if self.TeX_switch is False :
			leg_entries.append([histos['bgtot'], self.get_processName('bgtot'), 'l'])
		leg_entries.append([histos['pred' ], 'BG uncertainty'           , 'fl'])
		leg = self.draw_legend(leg_entries)

		# set minimum and maximum
		scale = 1.8
		if var == 'Int' : scale = 2.4
		if var.endswith('TotalBin') : scale = 1.4
		maximum = self.get_maximum(histos.values() + [hs_pred,], scale = scale, set_maximum = True, minimum = 0.)
#		gr_obs .SetMaximum(maximum)
#		gr_obs .SetMinimum(0.)

		# data histogram settings
#		self.apply_histoStyle(histos['obs'  ], 0)
#		self.apply_histoStyle(gr_obs         , 0)

		# special settings for total BG line and uncertainty
		histos['bgtot'].SetLineWidth(3)
		histos['bgtot'].SetLineColor(1)
#		histos['bgtot'].SetFillColor(12)
		histos['bgtot'].SetFillStyle(0)

		histos['pred' ].SetLineWidth(0)
		histos['pred' ].SetMarkerSize(0)
#		histos['pred' ].SetFillColor(1)
		histos['pred' ].SetFillColor(12)
		histos['pred' ].SetFillStyle(3005)

		# draw and set axis titles
		hs_pred.Draw('hist')
		self.set_axisTitles(hs_pred, self.get_varName(var), 'Events')
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
		elif var.startswith('CFChan') :
			total_bin = 0
			if var == 'CFChan_TotalBin' :
				total_bin = 1
				hs_pred.GetXaxis().SetBinLabel(total_bin, 'Total')
			for ibin in range(1 + total_bin, hs_pred.GetXaxis().GetNbins()+1) :
				bin = ibin - total_bin
				if self.TeX_switch is True :
					if bin < 4 : charge_str = 'p'
					else       : charge_str = 'm'
					if bin == 1 or bin == 4 : binlabel = '\\PGm'+charge_str+'\\PGm'+charge_str
					if bin == 2 or bin == 5 : binlabel = '\\Pe'+charge_str+'\\PGm'+charge_str
					if bin == 3 or bin == 6 : binlabel = '\\Pe'+charge_str+'\\Pe'+charge_str
				else :
					if bin < 4 : charge_str = '^{+}'
					else       : charge_str = '^{-}'
					binlabel = '#color[0]{l}'
					if bin == 1 or bin == 4 : binlabel += '#mu'+charge_str+'#mu'+charge_str
					if bin == 2 or bin == 5 : binlabel += 'e'+charge_str+'#mu'+charge_str
					if bin == 3 or bin == 6 : binlabel += 'e'+charge_str+'e'+charge_str
					binlabel += '#color[0]{l}'
				hs_pred.GetXaxis().SetBinLabel(ibin, binlabel)
			hs_pred.GetXaxis().SetLabelSize(0.062)
		elif var == 'Charge' :
			for bin in range(1, hs_pred.GetXaxis().GetNbins()+1) :
				if bin == 1 : binlabel = '--'
				if bin == 2 : binlabel = '++'
				hs_pred.GetXaxis().SetBinLabel(bin, binlabel)
				hs_pred.GetXaxis().SetLabelSize(0.062)
		elif not hs_pred.GetXaxis().IsVariableBinSize() and var != 'NVrtx' :
			bin_width = hs_pred.GetXaxis().GetBinWidth(1)
			hs_pred.GetYaxis().SetTitle(self.get_eventsPerGeVString(bin_width))
		leg.Draw()
		histos['pred'].Draw('0 E2 same')
		if self.TeX_switch is False :
			histos['bgtot'].Draw('hist same')
		if self.asymmErr :
#			gr_obs         .Draw('PE same')
			histos['obs'  ].SetBinErrorOption(ROOT.TH1.kPoisson)
			histos['obs'  ].Draw('PE0 X0 same')
		else :
			histos['obs'  ].Draw('PE X0 same')
		self.draw_cmsLine()
#		raw_input('ok? ')

		canvas.Update()
		self.save_canvas(canvas, self.path, 'ObsPred%s_%s%s' % (prefix, var, suffix))


	def save_plot_1d(self, h_data, h_mc, name = '', x_title = '', y_title = '') :
		if name == '' : name = h_mc.GetName()
		h_data = h_data.Clone()
		h_mc   = h_mc  .Clone()
		self.apply_histoStyle(h_data, 0)
		self.apply_histoStyle(h_mc  , 1)
		scale = 1.4
		maximum = scale * max(h_data.GetMaximum(), h_mc.GetMaximum())
		h_data.SetMaximum(maximum)
		h_mc  .SetMaximum(maximum)
		canvas = self.get_canvas()
		canvas.cd()
		leg_entries = [[h_data, 'Data', 'lp'], [h_mc, 'Simulation', 'lp']]
		leg = self.draw_legend(leg_entries)
		h_mc.Draw()
		self.set_axisTitles(h_mc, x_title, y_title)
		h_data.Draw('same pe')
		self.draw_cmsLine()
		leg.Draw()
		self.save_canvas(canvas, self.path, name)


	def save_plot_2d(self, histo, name = '', x_title = '', y_title = '') :
		if name == '' : name = histo.GetName()
		canvas = self.get_canvas()
		canvas.cd()
		histo.Draw('colztext')
		self.set_axisTitles(histo, x_title, y_title)
		self.draw_cmsLine()
		canvas.Print('%s%s.pdf' % (self.path, name))
		canvas.Print('%s%s.png' % (self.path, name))


#	def get_maximum(self, histos, scale = 1.8, set_maximum = True) :
#		'''returns maximum of a list of histograms'''
#		maximum = scale * max([histo.GetMaximum() for histo in histos])
#		if set_maximum :
#			for histo in histos :
#				histo.SetMaximum(maximum)
#		return maximum


	def drawTopLine(self, rightedge = 0.60, scale = 1., leftedge = 0.13) :
		latex = ROOT.TLatex()
		latex.SetNDC()
		latex.SetTextFont(62)
		latex.SetTextSize(scale * 0.05)
		latex.DrawLatex(leftedge, 0.92, 'CMS Preliminary')
		latex.SetTextFont(42)
		latex.SetTextSize(scale * 0.04)
		if self.lumi > 500. :
			lumi = self.lumi/1000.
			unit = 'fb^{-1}'
		else :
			lumi = self.lumi
			unit = 'pb^{-1}'
		latex.DrawLatex(rightedge, 0.92, '#sqrt{s} = 8 TeV, L_{int} = %4.1f %s' % (lumi, unit))


	def apply_histoStyle(self, histo, datamc) :
		if   datamc == 0 or datamc == 'data' or datamc == 'obs' :
			histo.SetMarkerStyle(20)
			histo.SetLineColor(ROOT.kBlack)
		elif datamc == 1 :
			histo.SetMarkerStyle(23)
			histo.SetMarkerColor(ROOT.kRed)
			histo.SetLineColor(ROOT.kRed)
		elif datamc == 'wjets' :
			histo.SetFillColor(ROOT.kOrange)
		elif datamc == 'zjets' :
			histo.SetFillColor(ROOT.kGreen)
		elif datamc == 'qcd' :
			histo.SetFillColor(ROOT.kYellow)
		else :
			histo.SetFillStyle(0)
		histo.SetMarkerSize(1.1)
		histo.SetLineWidth(2)
		histo.SetMinimum(0.)


	def set_axisTitles(self, histo, x_axis_title, y_axis_title) :
		histo.GetXaxis().SetTitle(x_axis_title)
		histo.GetYaxis().SetTitle(y_axis_title)
#		histo.GetXaxis().SetTitleOffset(1.25)
#		histo.GetYaxis().SetTitleOffset(1.25)
#		histo.GetXaxis().SetTitleSize(0.046)
#		histo.GetYaxis().SetTitleSize(0.046)
#		histo.GetXaxis().SetLabelSize(0.04)
#		histo.GetYaxis().SetLabelSize(0.04)
#		histo.GetXaxis().SetLabelOffset(0.012)
#		histo.GetYaxis().SetLabelOffset(0.012)


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
