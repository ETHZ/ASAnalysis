#! /usr/bin/python
import ROOT
import os, sys
import helper
import ttvStyle


class ttvplot(object) :

	def __init__(self, path, chan, lumi = 19500., cms_label = 0, asymmErr = True, TeX_switch = False, short_names = False) :
		self.ttvStyle = ttvStyle.ttvStyle(lumi = lumi, cms_label = cms_label, TeX_switch = TeX_switch)
		self.path = path
		if not self.path.endswith('/') : self.path += '/'
		if not os.path.exists(self.path) :
			os.makedirs(self.path)
		self.chan = chan
		self.asymmErr = asymmErr

		# random variable
		self.rand = ROOT.TRandom3(0)


	def save_plot(self, histos, var, prefix = '', suffix = '', charge_str = '') :
		'''save plot with observation and predictions'''

#		if selname != '' : selname += '_'
		if prefix  != '' and not prefix.startswith('_') : prefix = '_' + prefix
		if suffix  != '' and not suffix.startswith('_') : suffix = '_' + suffix

		# data with asymmetric errors
		gr_obs = helper.getGraphPoissonErrors_new(histos['obs'])

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
		if var.endswith('TotalBin') : scale = 1.4
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
			histo.SetFillColor(self.ttvStyle.get_fillColor(process))

		# data histogram settings
		self.apply_histoStyle(histos['obs'  ], 0)
		self.apply_histoStyle(gr_obs         , 0)
		
		# special settings for total BG line and uncertainty
		histos['bgtot'].SetLineWidth(3)
		histos['bgtot'].SetLineColor(1)
#		histos['bgtot'].SetFillColor(12)
		histos['bgtot'].SetFillStyle(0)

		histos['pred' ].SetLineWidth(0)
		histos['pred' ].SetFillColor(12)
		histos['pred' ].SetFillStyle(3005)

		# legend
		leg_entries = []
		leg_entries.append([histos['obs'  ], self.ttvStyle.get_processName('obs'  ), 'lp'])
		if self.chan == '3L' :
			leg_entries.append([histos['ttz'  ], self.ttvStyle.get_processName('ttz'  ), 'f'])
			leg_entries.append([histos['ttw'  ], self.ttvStyle.get_processName('ttw'  ), 'f'])
			leg_entries.append([histos['rare' ], self.ttvStyle.get_processName('rare' ), 'f'])
			leg_entries.append([histos['btag' ], self.ttvStyle.get_processName('btag' ), 'f'])
		if self.chan == '4L' :
			leg_entries.append([histos['ttz'  ], self.ttvStyle.get_processName('ttz'  ), 'f'])
			leg_entries.append([histos['ttw'  ], self.ttvStyle.get_processName('ttw'  ), 'f'])
			leg_entries.append([histos['rare' ], self.ttvStyle.get_processName('rare' ), 'f'])
			leg_entries.append([histos['zz'   ], self.ttvStyle.get_processName('zz'   ), 'f'])
		if self.chan == '2L' :
			leg_entries.append([histos['ttw'  ], self.ttvStyle.get_processName('ttw'  ), 'f'])
			leg_entries.append([histos['ttz'  ], self.ttvStyle.get_processName('ttz'  ), 'f'])
			leg_entries.append([histos['wz'   ], self.ttvStyle.get_processName('wz'   ), 'f'])
			leg_entries.append([histos['rare' ], self.ttvStyle.get_processName('rare' ), 'f'])
			leg_entries.append([histos['chmid'], self.ttvStyle.get_processName('chmid'), 'f'])
		leg_entries.append([histos['fake' ], self.ttvStyle.get_processName('fake' ), 'f'])
		leg_entries.append([histos['bgtot'], self.ttvStyle.get_processName('bgtot'), 'l'])
		leg_entries.append([histos['pred' ], 'BG uncertainty'           , 'fl'])
		leg = self.ttvStyle.draw_legend(leg_entries)

		canvas = self.ttvStyle.get_canvas('C_ObsPred')
		canvas.cd()

#		ROOT.gStyle.SetOptStat(0)
		ROOT.gStyle.SetOptTitle(0)
		ROOT.gStyle.SetEndErrorSize(0)  # set the size of the small line at the end of the error bars
		ROOT.gPad.SetTicks(1,1)

		# draw and set axis titles
		hs_pred.Draw()
		self.set_axisTitles(hs_pred, self.ttvStyle.get_varName(var), 'Events')
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
			hs_pred.GetYaxis().SetTitle('Events / %.0f GeV' % bin_width)
		leg.Draw()
		histos['pred'].Draw('0 E2 same')
		histos['bgtot'].Draw('hist same')
		if self.asymmErr :
			gr_obs         .Draw('PE same')
		else :
			histos['obs'  ].Draw('PE X0 same')
		self.ttvStyle.draw_cmsLine()
#		raw_input('ok? ')

		canvas.Print('%sObsPred%s_%s%s.pdf' % (self.path, prefix, var, suffix))
		canvas.Print('%sObsPred%s_%s%s.png' % (self.path, prefix, var, suffix))
#		raw_input('ok? ')


	def save_controlPlot(self, histos, var, prefix = '') :
		if prefix  != '' : prefix += '_'
		canvas = self.ttvStyle.get_canvas()
		canvas.cd()
		hstack = ROOT.THStack('hs_%s' % var, '%s' % var)
		leg_entries = []
		for process in histos :
			self.apply_histoStyle(histos[process], process)
			if process == 'data' or process == 'obs' :
				leg_entries.append([histos[process], self.ttvStyle.get_processName('obs'), 'lp'])
				data_index = len(leg_entries) - 1
			else :
				hstack.Add(histos[process])
				leg_entries.append([histos[process], self.ttvStyle.get_processName(process), 'f'])
		leg_entries[0], leg_entries[data_index] = leg_entries[data_index], leg_entries[0]
		maximum = self.get_maximum(histos.values())
		hstack.SetMaximum(maximum)
		leg = self.ttvStyle.draw_legend(leg_entries)
		hstack.Draw('hist')
		bin_width = hstack.GetXaxis().GetBinWidth(1)
		y_title = 'Events / %.0f GeV' % bin_width
		self.set_axisTitles(hstack, self.ttvStyle.get_varName(var), y_title)
		histos['data'].Draw('psame')
		self.ttvStyle.draw_cmsLine()
		canvas.cd()
		leg.Draw()
		canvas.Print('%s%s%s.pdf' % (self.path, prefix, var))
		canvas.Print('%s%s%s.png' % (self.path, prefix, var))


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
		canvas = self.ttvStyle.get_canvas()
		canvas.cd()
		leg_entries = [[h_data, 'Data', 'lp'], [h_mc, 'Simulation', 'lp']]
		leg = self.ttvStyle.draw_legend(leg_entries)
		h_mc.Draw()
		self.set_axisTitles(h_mc, x_title, y_title)
		h_data.Draw('same pe')
		self.ttvStyle.draw_cmsLine()
		leg.Draw()
		canvas.Print('%s%s.pdf' % (self.path, name))
		canvas.Print('%s%s.png' % (self.path, name))


	def save_plot_2d(self, histo, name = '', x_title = '', y_title = '') :
		if name == '' : name = histo.GetName()
		canvas = self.ttvStyle.get_canvas()
		canvas.cd()
		histo.Draw('colztext')
		self.set_axisTitles(histo, x_title, y_title)
		self.ttvStyle.draw_cmsLine()
		canvas.Print('%s%s.pdf' % (self.path, name))
		canvas.Print('%s%s.png' % (self.path, name))


	def get_maximum(self, histos, scale = 1.8, set_maximum = True) :
		'''returns maximum of a list of histograms'''
		maximum = scale * max([histo.GetMaximum() for histo in histos])
		if set_maximum :
			for histo in histos :
				histo.SetMaximum(maximum)
		return maximum


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
