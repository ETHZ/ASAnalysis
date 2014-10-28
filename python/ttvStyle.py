#! /usr/bin/python
import ROOT
from array import array


class ttvStyle(object) :

	def __init__(self, lumi = 19500., cms_label = 0, TeX_switch = False, short_names = False) :
		self.ttvStyle = ROOT.TStyle('ttvStyle','ttV Style')
		self.lumi = lumi
		self.cms_label = cms_label
		self.TeX_switch = TeX_switch
		self.short_names = short_names

		# random variable
		self.rand = ROOT.TRandom3(0)

		self.set_style()
		self.define_colors()
		self.define_processNames(self.TeX_switch)
		self.define_varNames()


	def set_style(self) :
		'''make all the style settings'''

		# canvas
		self.ttvStyle.SetCanvasBorderMode(0)
		self.ttvStyle.SetCanvasColor(ROOT.kWhite)
		self.ttvStyle.SetCanvasDefH(600)
		self.ttvStyle.SetCanvasDefW(600)
		self.ttvStyle.SetCanvasDefX(0)
		self.ttvStyle.SetCanvasDefY(0)

		# pad
#		self.ttvStyle.SetPadBorderMode(0)
#		# self.ttvStyle.SetPadBorderSize(Width_t size = 1)
#		self.ttvStyle.SetPadColor(ROOT.kWhite)
#		self.ttvStyle.SetPadGridX(False)
#		self.ttvStyle.SetPadGridY(False)
#		self.ttvStyle.SetGridColor(0)
#		self.ttvStyle.SetGridStyle(3)
#		self.ttvStyle.SetGridWidth(1)

		# frame
		self.ttvStyle.SetFrameBorderMode(0)
		self.ttvStyle.SetFrameBorderSize(1)
		self.ttvStyle.SetFrameFillColor(0)
		self.ttvStyle.SetFrameFillStyle(0)
		self.ttvStyle.SetFrameLineColor(1)
		self.ttvStyle.SetFrameLineStyle(1)
		self.ttvStyle.SetFrameLineWidth(1)

		# histo
		self.ttvStyle.SetHistFillColor(63)
		self.ttvStyle.SetHistFillStyle(1001) # 0: hollow, 1001: solid
		self.ttvStyle.SetHistLineColor(1)
		self.ttvStyle.SetHistLineStyle(0)
		self.ttvStyle.SetHistLineWidth(1)
#		# self.ttvStyle.SetLegoInnerR(Float_t rad = 0.5)
#		# self.ttvStyle.SetNumberContours(Int_t number = 20)
#
#		self.ttvStyle.SetEndErrorSize(2)
##		self.ttvStyle.SetErrorMarker(20)
#		#self.ttvStyle.SetErrorX(0.)
#	
		self.ttvStyle.SetMarkerStyle(20)
		self.ttvStyle.SetMarkerColor(1)
#		self.ttvStyle.SetMarkerSize(1.2)
#
#		# fit/function
##		self.ttvStyle.SetOptFit(1)
		self.ttvStyle.SetOptFit(0)
#		self.ttvStyle.SetFitFormat('5.4g')
		self.ttvStyle.SetFuncColor(3)
#		self.ttvStyle.SetFuncStyle(1)
#		self.ttvStyle.SetFuncWidth(1)
#
#		# date
#		self.ttvStyle.SetOptDate(0)
#		# self.ttvStyle.SetDateX(Float_t x = 0.01)
#		# self.ttvStyle.SetDateY(Float_t y = 0.01)

		# statistics box
#		self.ttvStyle.SetOptFile(0)
		self.ttvStyle.SetOptStat(0)  # To display the mean and RMS:   SetOptStat('mr')
		self.ttvStyle.SetStatColor(ROOT.kWhite)
		self.ttvStyle.SetStatFont(42)
		self.ttvStyle.SetStatFontSize(0.025)
		self.ttvStyle.SetStatTextColor(1)
#		self.ttvStyle.SetStatFormat('6.4g')
		self.ttvStyle.SetStatBorderSize(0)
#		self.ttvStyle.SetStatH(0.1)
		self.ttvStyle.SetStatW(0.35)
		self.ttvStyle.SetStatStyle(1001)
		self.ttvStyle.SetStatX(0.9)
		self.ttvStyle.SetStatY(0.9)

		# margins
		self.ttvStyle.SetPadTopMargin   (0.06)
		self.ttvStyle.SetPadBottomMargin(0.14)
		self.ttvStyle.SetPadLeftMargin  (0.14)
		self.ttvStyle.SetPadRightMargin (0.06)

		# global title
		self.ttvStyle.SetOptTitle(0)
		self.ttvStyle.SetTitleFont(42)
		self.ttvStyle.SetTitleColor(1)
#		self.ttvStyle.SetTitleTextColor(1)
		self.ttvStyle.SetTitleFillColor(10)
#		self.ttvStyle.SetTitleFontSize(0.05)
#		# self.ttvStyle.SetTitleH(0)  # Set the height of the title box
#		# self.ttvStyle.SetTitleW(0)  # Set the width of the title box
#		# self.ttvStyle.SetTitleX(0)  # Set the position of the title box
#		# self.ttvStyle.SetTitleY(0.985)  # Set the position of the title box
#		# self.ttvStyle.SetTitleStyle(Style_t style = 1001)
#		# self.ttvStyle.SetTitleBorderSize(2)
#
#		# axis titles
#		self.ttvStyle.SetTitleColor(1, 'XYZ')
		self.ttvStyle.SetTitleFont(42, 'XYZ')
		self.ttvStyle.SetTitleSize(0.046, 'XYZ')
#		self.ttvStyle.SetTitleXSize(0.02)  # Another way to set the size?
#		self.ttvStyle.SetTitleYSize(0.02)
#		self.ttvStyle.SetTitleXOffset(1.25)
#		self.ttvStyle.SetTitleYOffset(1.25)
		self.ttvStyle.SetTitleOffset(1.5, 'XYZ')  # Another way to set the Offset
#
#		# axis labels
#		self.ttvStyle.SetLabelColor(1, 'XYZ')
		self.ttvStyle.SetLabelFont(42, 'XYZ')
		self.ttvStyle.SetLabelOffset(0.012, 'XYZ')
		self.ttvStyle.SetLabelSize(0.04, 'XYZ')
#
#		# axis
#		self.ttvStyle.SetAxisColor(1, 'XYZ')
#		self.ttvStyle.SetStripDecimals(ROOT.kTRUE)
#		self.ttvStyle.SetTickLength(0.03, 'XYZ')
#		self.ttvStyle.SetNdivisions(510, 'XYZ')
		self.ttvStyle.SetPadTickX(1)   # To get tick marks on the opposite side of the frame
		self.ttvStyle.SetPadTickY(1)
#
#		# Change for log plots:
#		self.ttvStyle.SetOptLogx(0)
#		self.ttvStyle.SetOptLogy(0)
#		self.ttvStyle.SetOptLogz(0)
#
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
		if TeX_switch is True :
			self.process_names['ttz'  ] = '\\ttz'
			self.process_names['ttw'  ] = '\\ttw'
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
		self.var_names['pT1'   ] = 'Leading lepton p_{T} [GeV]'
		self.var_names['pT2'   ] = 'Sub-leading lepton p_{T} [GeV]'
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


	def get_fillColor(self, process) :
		if process in self.colors.keys() :
			return self.colors[process]
		return 2


	def get_varName(self, var) :
		if var in self.var_names.keys() :
			return self.var_names[var]
		return ''


	def get_processName(self, process) :
		if process in self.process_names.keys() :
			return self.process_names[process]
		return '?'


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


	def draw_legend(self, entries) :
		leg = ROOT.TLegend()
		for entry in entries :
			leg.AddEntry(entry[0], ' ' + entry[1], entry[2])

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
		if self.lumi > 500. :
			lumi = self.lumi/1000.
			unit = 'fb^{-1}'
		else :
			lumi = self.lumi
			unit = 'pb^{-1}'
#		latex.DrawLatex(0.15, 0.88, '%4.1f %s (8 TeV)' % (lumi, unit))
		latex.SetTextAlign(31)
		latex.DrawLatex(lumi_x, lumi_y, '%4.1f %s (8 TeV)' % (lumi, unit))


	def get_maximum(self, histos, scale = 1.8, set_maximum = True) :
		'''returns maximum of a list of histograms'''
		maximum = scale * max([histo.GetMaximum() for histo in histos])
		if set_maximum :
			for histo in histos :
				histo.SetMaximum(maximum)
		return maximum
