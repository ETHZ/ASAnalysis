#! /usr/bin/python
import ROOT
from array import array


def ttvStyle() :
	ttvStyle = ROOT.TStyle('ttvStyle','RD42 Style')

	# canvas
	ttvStyle.SetCanvasBorderMode(0)
	ttvStyle.SetCanvasColor(ROOT.kWhite)
	ttvStyle.SetCanvasDefH(600)
	ttvStyle.SetCanvasDefW(600)
	ttvStyle.SetCanvasDefX(0)
	ttvStyle.SetCanvasDefY(0)

	# pad
#	ttvStyle.SetPadBorderMode(0)
#	# ttvStyle.SetPadBorderSize(Width_t size = 1)
#	ttvStyle.SetPadColor(ROOT.kWhite)
#	ttvStyle.SetPadGridX(False)
#	ttvStyle.SetPadGridY(False)
#	ttvStyle.SetGridColor(0)
#	ttvStyle.SetGridStyle(3)
#	ttvStyle.SetGridWidth(1)

	# frame
	ttvStyle.SetFrameBorderMode(0)
	ttvStyle.SetFrameBorderSize(1)
	ttvStyle.SetFrameFillColor(0)
	ttvStyle.SetFrameFillStyle(0)
	ttvStyle.SetFrameLineColor(1)
	ttvStyle.SetFrameLineStyle(1)
	ttvStyle.SetFrameLineWidth(1)

	# histo
##	ttvStyle.SetHistFillColor(63)
#	# ttvStyle.SetHistFillStyle(0)
#	ttvStyle.SetHistLineColor(1)
#	ttvStyle.SetHistLineStyle(0)
#	ttvStyle.SetHistLineWidth(1)
#	# ttvStyle.SetLegoInnerR(Float_t rad = 0.5)
#	# ttvStyle.SetNumberContours(Int_t number = 20)
#
#	ttvStyle.SetEndErrorSize(2)
##	ttvStyle.SetErrorMarker(20)
#	#ttvStyle.SetErrorX(0.)
#	
#	ttvStyle.SetMarkerStyle(20)
#	ttvStyle.SetMarkerSize(1.2)
#
#	# fit/function
##	ttvStyle.SetOptFit(1)
	ttvStyle.SetOptFit(0)
#	ttvStyle.SetFitFormat('5.4g')
	ttvStyle.SetFuncColor(3)
#	ttvStyle.SetFuncStyle(1)
#	ttvStyle.SetFuncWidth(1)
#
#	# date
#	ttvStyle.SetOptDate(0)
#	# ttvStyle.SetDateX(Float_t x = 0.01)
#	# ttvStyle.SetDateY(Float_t y = 0.01)

	# statistics box
#	ttvStyle.SetOptFile(0)
	ttvStyle.SetOptStat(0)  # To display the mean and RMS:   SetOptStat('mr')
	ttvStyle.SetStatColor(ROOT.kWhite)
	ttvStyle.SetStatFont(42)
	ttvStyle.SetStatFontSize(0.025)
	ttvStyle.SetStatTextColor(1)
#	ttvStyle.SetStatFormat('6.4g')
	ttvStyle.SetStatBorderSize(0)
#	ttvStyle.SetStatH(0.1)
	ttvStyle.SetStatW(0.35)
	ttvStyle.SetStatStyle(1001)
	ttvStyle.SetStatX(0.9)
	ttvStyle.SetStatY(0.9)

	# margins
	ttvStyle.SetPadTopMargin   (0.06)
	ttvStyle.SetPadBottomMargin(0.14)
	ttvStyle.SetPadLeftMargin  (0.14)
	ttvStyle.SetPadRightMargin (0.06)

	# global title
	ttvStyle.SetOptTitle(0)
	ttvStyle.SetTitleFont(42)
	ttvStyle.SetTitleColor(1)
#	ttvStyle.SetTitleTextColor(1)
	ttvStyle.SetTitleFillColor(10)
#	ttvStyle.SetTitleFontSize(0.05)
#	# ttvStyle.SetTitleH(0)  # Set the height of the title box
#	# ttvStyle.SetTitleW(0)  # Set the width of the title box
#	# ttvStyle.SetTitleX(0)  # Set the position of the title box
#	# ttvStyle.SetTitleY(0.985)  # Set the position of the title box
#	# ttvStyle.SetTitleStyle(Style_t style = 1001)
#	# ttvStyle.SetTitleBorderSize(2)
#
#	# axis titles
#	ttvStyle.SetTitleColor(1, 'XYZ')
	ttvStyle.SetTitleFont(42, 'XYZ')
	ttvStyle.SetTitleSize(0.046, 'XYZ')
#	ttvStyle.SetTitleXSize(0.02)  # Another way to set the size?
#	ttvStyle.SetTitleYSize(0.02)
#	ttvStyle.SetTitleXOffset(1.25)
#	ttvStyle.SetTitleYOffset(1.25)
	ttvStyle.SetTitleOffset(1.5, 'XYZ')  # Another way to set the Offset
#
#	# axis labels
#	ttvStyle.SetLabelColor(1, 'XYZ')
	ttvStyle.SetLabelFont(42, 'XYZ')
	ttvStyle.SetLabelOffset(0.012, 'XYZ')
	ttvStyle.SetLabelSize(0.04, 'XYZ')
#
#	# axis
#	ttvStyle.SetAxisColor(1, 'XYZ')
#	ttvStyle.SetStripDecimals(ROOT.kTRUE)
#	ttvStyle.SetTickLength(0.03, 'XYZ')
#	ttvStyle.SetNdivisions(510, 'XYZ')
	ttvStyle.SetPadTickX(1)   # To get tick marks on the opposite side of the frame
	ttvStyle.SetPadTickY(1)
#
#	# Change for log plots:
#	ttvStyle.SetOptLogx(0)
#	ttvStyle.SetOptLogy(0)
#	ttvStyle.SetOptLogz(0)
#
#	# colors
	ttvStyle.SetPalette(1)
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
	ttvStyle.SetNumberContours(ncontours)

	# text
#	ttvStyle.SetTextAlign(12)

#
#	# postscript options:
#	#ttvStyle.SetPaperSize(20.,20.)
##	ttvStyle.SetLineScalePS(Float_t scale = 3)
##	ttvStyle.SetLineStyleString(Int_t i, const char* text)
##	ttvStyle.SetHeaderPS(const char* header)
##	ttvStyle.SetTitlePS(const char* pstitle)
#
##	ttvStyle.SetBarOffset(Float_t baroff = 0.5)
##	ttvStyle.SetBarWidth(Float_t barwidth = 0.5)
##	ttvStyle.SetPaintTextFormat(const char* format = 'g')
##	ttvStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
##	ttvStyle.SetTimeOffset(Double_t toffset)
##	ttvStyle.SetHistMinimumZero(kTRUE)

	ttvStyle.cd()
