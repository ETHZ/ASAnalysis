#! /usr/bin/python
import ROOT
from array import array
import helper
import ConfigParser
import os
import sys


class ttvStyle(object) :

	def __init__(self, lumi = -1., cms_label = 0, TeX_switch = False, short_names = False, config_file = '') :
		if config_file == '' :
			config_file = '%s/ttvStyle.cfg' % os.path.dirname(os.path.realpath(__file__))
		self.ttvStyle = self.parse_styleSettings(config_file, style_name = 'ttvStyle', style_title = 'ttV Style')
		self.lumi = lumi
		self.cms_label = cms_label
		self._TeX_switch = TeX_switch
		self.short_names = short_names

		# random variable
		self.rand = ROOT.TRandom3(0)

		self.set_style()
		self.define_colors()
		self.define_processNames(self.TeX_switch)
		self.define_varNames()


	@property
	def TeX_switch(self) :
		return self._TeX_switch


	def parse_styleSettings(self, config_file, style_name = 'ttvStyle', style_title = 'ttV Style') :
		if not os.path.exists(config_file) :
			print '[ERROR] %s does not exist!' % config_file
			sys.exit(1)
		style = ROOT.TStyle(style_name, style_title)
		config = ConfigParser.ConfigParser()
		config.optionxform = str # case sensitive options
		config.read(config_file)
		for section in config.sections() :
			for key, value_str in config.items(section) :
				value = eval(value_str)
				if type(value) is not tuple :
					getattr(style, 'Set%s' % key)(value)
				else :
					getattr(style, 'Set%s' % key)(*value)
		return style


	def set_style(self) :
		'''make all the style settings'''

		# PostScript output
		if self.TeX_switch : self.ttvStyle.SetPaperSize(15., 15.)

#		# colors
		self.ttvStyle.SetPalette(1)
		stops = [0.00, 0.34, 0.61, 0.84, 1.00]
		red   = [0.00, 0.00, 0.87, 1.00, 0.51]
		green = [0.00, 0.81, 1.00, 0.20, 0.00]
		blue  = [0.51, 1.00, 0.12, 0.00, 0.00]
		s = array('d', stops)
		r = array('d', red)
		g = array('d', green)
		b = array('d', blue)
		ncontours = 999
		npoints = len(s)
		ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
		self.ttvStyle.SetNumberContours(ncontours)

		# text
#		self.ttvStyle.SetTextAlign(12)

#
#		# postscript options:
#		#self.ttvStyle.SetPaperSize(20.,20.)
##		self.ttvStyle.SetLineScalePS(Float_t scale = 3)
##		self.ttvStyle.SetLineStyleString(Int_t i, const char* text)
##		self.ttvStyle.SetHeaderPS(const char* header)
##		self.ttvStyle.SetTitlePS(const char* pstitle)
#
##		self.ttvStyle.SetBarOffset(Float_t baroff = 0.5)
##		self.ttvStyle.SetBarWidth(Float_t barwidth = 0.5)
##		self.ttvStyle.SetPaintTextFormat(const char* format = 'g')
##		self.ttvStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
##		self.ttvStyle.SetTimeOffset(Double_t toffset)
##		self.ttvStyle.SetHistMinimumZero(kTRUE)

		self.ttvStyle.cd()


	def define_colors(self) :
		'''define histogram colors'''

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
		self.colors['wjets'] = 36
		self.colors['zjets'] = 49
		self.colors['qcd'  ] = 40
		self.colors['top'  ] = 46


	def define_processNames(self, TeX_switch = False) :
		'''define process names'''

		self.process_names = {}
		self.process_names['obs'  ] = 'Observed'
		self.process_names['fake' ] = 'Misidentified lepton'
		self.process_names['chmid'] = 'Mismeasured charge'
		self.process_names['rare' ] = 'Irreducible'
		self.process_names['wz'   ] = 'WZ'
		self.process_names['ttz'  ] = 't#bar{t}Z'
		self.process_names['ttw'  ] = 't#bar{t}W'
		self.process_names['bgtot'] = 'Backgrounds'
		self.process_names['btag' ] = 'Non-top'
		self.process_names['zz'   ] = 'ZZ'
		self.process_names['wjets'] = 'W+jets'
		self.process_names['zjets'] = 'DY+jets'
		self.process_names['qcd'  ] = 'QCD'
		self.process_names['top'  ] = 'Top'
		self.process_names['ttbar'] = 't#bar{t}'
		if TeX_switch is True :
			self.process_names['obs'  ] = '\\data'
			self.process_names['fake' ] = '\\fakes'
			self.process_names['chmid'] = '\\cmid'
			self.process_names['rare' ] = '\\irreducible'
			self.process_names['ttz'  ] = '\\ttz'
			self.process_names['ttw'  ] = '\\ttw'
			self.process_names['wz'   ] = '\\wz'
			self.process_names['wjets'] = '\\wjets'
			self.process_names['zjets'] = '\\zjets'
			self.process_names['qcd'  ] = '\\qcd'
			self.process_names['ttbar'] = '\\ttbar'
			if self.short_names is True :
				self.process_names['fake' ] = 'Lept.\\ MisID'
				self.process_names['chmid'] = 'Ch.\\ MisID'


	def define_varNames(self) :
		''' define variable names'''

		self.var_names = {}
		self.var_names['HT'    ] = 'H_{T} [GeV]'
		self.var_names['MET'   ] = 'Particle flow E_{T}^{miss} [GeV]'
		self.var_names['NJ'    ] = 'Jet multiplicity'
		self.var_names['NbJ'   ] = 'b-Jet multiplicity (CSVL)'
		self.var_names['NbJmed'] = 'b-Jet multiplicity (CSVM)'
		self.var_names['pT'    ] = 'p_{T} [GeV]'
		self.var_names['pT1'   ] = 'Leading lepton p_{T} [GeV]'
		self.var_names['pT2'   ] = 'Sub-leading lepton p_{T} [GeV]'
		self.var_names['pTl'   ] = 'Lepton p_{T} [GeV]'
		self.var_names['pTm'   ] = 'Muon p_{T} [GeV]'
		self.var_names['pTe'   ] = 'Electron p_{T} [GeV]'
		self.var_names['eta'   ] = '#eta'
		self.var_names['eta1'  ] = 'Leading lepton #eta'
		self.var_names['eta2'  ] = 'Sub-leading lepton #eta'
		self.var_names['etal'  ] = 'Lepton #eta'
		self.var_names['etam'  ] = 'Muon #eta'
		self.var_names['etae'  ] = 'Electron #eta'
		self.var_names['Int'   ] = ''
		self.var_names['Mll'   ] = 'm_{ll} [GeV]'
		self.var_names['NVrtx' ] = 'N_{Vertices}'
		self.var_names['minMT' ] = 'M_{T} [GeV]'
		self.var_names['M3'    ] = 'm_{bjj} [GeV]'
		self.var_names['NJets'      ] = self.get_varName('NJ')
		self.var_names['MaxJPt'     ] = 'Hardest jet p_{T} [GeV]'
		self.var_names['NVertices'  ] = self.get_varName('NVrtx')
		self.var_names['ClosJetPt'  ] = 'Closest jet p_{T} [GeV]'
		self.var_names['AwayJetPt'  ] = 'Away jet p_{T} [GeV]'
		self.var_names['NBJets'     ] = self.get_varName('NbJ')
		self.var_names['MT'         ] = 'm_{T}'
		self.var_names['MET_noMTCut'] = 'E_{T}^{miss} [GeV]'
		self.var_names['MT_MET30'   ] = self.get_varName('MT')
		self.var_names['LepPt'      ] = 'Lepton p_{T} [GeV]'
		self.var_names['LepEta'     ] = 'Lepton #eta [GeV]'
		self.var_names['LepIso'     ] = 'Lepton isolation'
		self.var_names['ClosJetDR'  ] = 'Closest jet DR'
		self.var_names['AwayJetDR'  ] = 'Away jet DR'
		self.var_names['fratio'     ] = 'f = N_{t}/N_{l}'
		self.var_names['pratio'     ] = 'p = N_{t}/N_{l}'
		self.var_names['effSig'     ] = '#varepsilon_{Signal} [%]'
		self.var_names['sigexp'     ] = '#sigma_{expected}'
		self.var_names['XSecErr'    ] = 'Exp. Cross Section Error [fb]'
		if self.TeX_switch is True :
			self.var_names['HT'    ] = '\\text{\\HT{} (\\si{\\giga\\electronvolt})}'
			self.var_names['MET'   ] = '\\text{\\MET{} (\\si{\\giga\\electronvolt})}'
			self.var_names['NJ'    ] = '\\nj'
			self.var_names['NbJ'   ] = '\\nb'
			self.var_names['NbJmed'] = '\\nb'
			self.var_names['pT'    ] = '\\text{\\pt{} (\\si{\\GeV})}'
			self.var_names['pT1'   ] = '\\text{Leading lepton \\pt{} (\\si{\\giga\\electronvolt})}'
			self.var_names['pT2'   ] = '\\text{Sub-leading lepton \pt{} (\\si{\\giga\\electronvolt})}'
			self.var_names['pTl'   ] = '\\text{Lepton \\pt{} (\\si{\\GeV})}'
			self.var_names['pTm'   ] = '\\text{Muon \\pt{} (\\si{\\GeV})}'
			self.var_names['pTe'   ] = '\\text{Electron \\pt{} (\\si{\\GeV})}'
			self.var_names['eta'   ] = '\\eta'
			self.var_names['eta1'  ] = '\\text{Leading lepton \\ensuremath{\\eta}}'
			self.var_names['eta2'  ] = '\\text{Sub-leading lepton \\ensuremath{\\eta}}'
			self.var_names['etal'  ] = '\\text{Lepton \\ensuremath{\\eta}}'
			self.var_names['etam'  ] = '\\text{Muon \\ensuremath{\\eta}}'
			self.var_names['etae'  ] = '\\text{Electron \\ensuremath{\\eta}}'
			self.var_names['Int'   ] = ''
			self.var_names['Mll'   ] = '\\text{\\mll{} (\\si{\\giga\\electronvolt})}'
			self.var_names['NVrtx' ] = '\\nvrtx'
			self.var_names['minMT' ] = '\\text{\\mt{} (\\si{\\giga\\electronvolt})}'
			self.var_names['M3'    ] = '\\text{\\mbjj{} (\\si{\\giga\\electronvolt})}'
			self.var_names['fratio'] = '\\fratio = \\nt / \\nl'
			self.var_names['pratio'] = '\\pratio = \\nt / \\nl'
			self.var_names['effSig'     ] = '\\text{\\effSig{} (\\si{\\percent})}'
			self.var_names['sigexp'     ] = '\\sigexp'
			self.var_names['XSecErr'    ] = '\\text{\\xsecexperr{} (\\si{\\femto\\barn})}'


	def get_fillColor(self, process) :
		if process in self.colors.keys() :
			return self.colors[process]
		return 2


	def get_fillStyle(self, process) :
		if process == 'ttw' :
			return 3004
		elif process == 'ttbar' :
			return 3005
		else :
			return 1001


	def get_varName(self, var) :
		if var in self.var_names.keys() :
			return self.var_names[var]
		return ''


	def get_processName(self, process) :
		if process in self.process_names.keys() :
			return self.process_names[process]
		return '?'


	def get_processCommand(self, process) :
		if   process == 'WWTo2L2Nu' : command = '\\ww'
		elif process == 'WZTo3LNu'  : command = '\\wz'
		elif process == 'ZZTo4L'    : command = '\\zz'
		else : command = '\\%s' % process.replace('+', 'p').replace('-', 'm').replace('_', '')
		return command


	def get_eventsPerGeVString(self, bin_width) :
		if self.TeX_switch :
			return '\\text{Events / \\SI{%g}{\\GeV}}' % bin_width
		else :
			return 'Events / %g GeV' % bin_width


	def label_binWidth(self, var) :
		variables = []
		variables.append('HT'    )
		variables.append('MET'   )
		variables.append('pT'    )
		variables.append('pT1'   )
		variables.append('pT2'   )
		variables.append('pTl'   )
		variables.append('pTm'   )
		variables.append('pTe'   )
		variables.append('Mll'   )
		variables.append('minMT' )
		variables.append('M3'    )
		variables.append('MaxJPt'     )
		variables.append('ClosJetPt'  )
		variables.append('AwayJetPt'  )
		variables.append('MT'         )
		variables.append('MET_noMTCut')
		variables.append('MT_MET30'   )
		variables.append('LepPt'      )
		if var in variables : return True
		else                : return False


	@property
	def cms_label(self) :
		if   self._cms_label == 0 : return 'CMS'
		elif self._cms_label == 1 : return 'CMS Simulation'
		elif self._cms_label == 2 : return 'CMS Preliminary'
		elif self._cms_label == 3 : return 'CMS Simulation Preliminary'
		else                      : return ''


	@cms_label.setter
	def cms_label(self, cms_label) :
		self._cms_label = cms_label


	def get_canvas(self, name = '') :
		name += '_%d' % self.rand.Integer(10000)  # add random number to avoid same names
		canvas = ROOT.TCanvas(name, name)
#		canvas.SetLeftMargin(0.12)
#		canvas.SetRightMargin(0.04)
#		canvas.SetTopMargin(0.04)
#		canvas.SetBottomMargin(0.12)
		canvas.cd()

#		ROOT.gStyle.SetOptStat(0)
#		ROOT.gStyle.SetOptTitle(0)
#		ROOT.gStyle.SetEndErrorSize(0)  # set the size of the small line at the end of the error bars
##		ROOT.gStyle.SetErrorX(0)
#		ROOT.gPad.SetTicks(1,1)
		return canvas


	def save_canvas(self, canvas, path, file_name) :
		if not path.endswith('/') : path += '/'
		helper.mkdir(path)
		canvas.Update()
		if self.TeX_switch :
			canvas.Print('%s%s.tex' % (path, file_name))
		else :
			canvas.Print('%s%s.pdf' % (path, file_name))
			canvas.Print('%s%s.png' % (path, file_name))
			canvas.Print('%s%s.root' % (path, file_name))


	def draw_legend(self, entries) :
		leg = ROOT.TLegend()
		if self.TeX_switch : extra_space = ''
		else               : extra_space = ' '
		for entry in entries :
			leg.AddEntry(entry[0], extra_space + entry[1], entry[2])

		# set position
		width = 0.17
		x = 1. - self.ttvStyle.GetPadRightMargin() - 0.33
		y = 1. - self.ttvStyle.GetPadTopMargin() - 0.03
		leg.SetX1NDC(x)
		leg.SetX2NDC(x+width)
		leg.SetY1NDC(y-leg.GetNRows()*0.25*width)
		leg.SetY2NDC(y)
		leg.SetFillStyle(0)
		leg.SetTextFont(42)
		leg.SetTextSize(0.03)
		leg.SetBorderSize(0)
		leg.SetTextAlign(12)

		leg.Draw()
		return leg


	def draw_cmsLine(self) :
		cms_x  = self.ttvStyle.GetPadLeftMargin() + 0.03
		cms_y  = 1. - self.ttvStyle.GetPadTopMargin() - 0.03
		lumi_x = 1. - self.ttvStyle.GetPadRightMargin()
		lumi_y = 1. - self.ttvStyle.GetPadTopMargin() + 0.01
		latex = ROOT.TLatex()
		latex.SetNDC()
		latex.SetTextFont(62)
		latex.SetTextSize(0.04)
		latex.SetTextAlign(13)
		latex.DrawLatex(cms_x, cms_y, self.cms_label)
		latex.SetTextFont(42)
		latex.SetTextSize(0.03)
		if self.lumi < 0 : return
		latex.SetTextAlign(31)
		if self.lumi > 500. :
			lumi = self.lumi/1000.
			unit = 'fb^{-1}'
			if self.TeX_switch is True :
				unit = '\\femto'
		else :
			lumi = self.lumi
			unit = 'pb^{-1}'
			if self.TeX_switch is True :
				unit = '\\pico'
#		latex.DrawLatex(0.15, 0.88, '%4.1f %s (8 TeV)' % (lumi, unit))
		if self.TeX_switch is True :
			latex.DrawLatex(lumi_x, lumi_y, '\\text{\\SI{%4.1f}{\\per%s\\barn} (\\SI{8}{\\tera\\electronvolt})}' % (lumi, unit))
		else :
			latex.DrawLatex(lumi_x, lumi_y, '%4.1f %s (8 TeV)' % (lumi, unit))


	@staticmethod
	def get_maximum(histos, scale = 1.8, set_maximum = True, minimum = None) :
		'''returns maximum of a list of histograms'''
		maximum = scale * max([histo.GetMaximum() for histo in histos])
		for histo in histos :
			if set_maximum :
				histo.SetMaximum(maximum)
			if minimum != None :
				histo.SetMinimum(minimum)
		return maximum
