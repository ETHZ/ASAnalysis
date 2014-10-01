#!/usr/bin/python
import ROOT
import os, sys
#from ttvStyle import ttvStyle
from array import array
import imp
ttvStyle = imp.load_source('ttvStyle', '../python/ttvStyle.py')


def make_table(central, pdfsets, sample = '') :
	'''print latex table'''

	if sample != '' : sample += '_'
	with open('%sacceptances.tex' % sample, 'w') as file :
		file.write('\\sisetup{\n')
		file.write('\ttable-number-alignment = center,\n')
		file.write('\ttable-figures-integer  = 1,\n')
		file.write('\ttable-figures-decimal  = 2,\n')
		file.write('\ttable-figures-exponent = 2\n')
		file.write('}\n\n')
		file.write('\\begin{tabular}{l|SSS}\n')
		file.write('\t\\hline\\hline\n')
		file.write('\t{PDF set}   & {$A_0$} & {$\\Delta A^+$} & {$\\Delta A^-$} \\\\\n')
		file.write('\t\\hline\n')
		file.write('\t\\pdfCTEQ{}  & %.2e & & \\\\\n' % central)
		file.write('\t\\hline\n')
		for pdfset in pdfsets :
			file.write('\t%s & %.2e & %.0e & %.0e \\\\\n' % pdfset[1:])
		file.write('\t\\hline\\hline\n')
		file.write('\\end{tabular}')


def make_plot(central, pdfsets, sample = '') :
	if sample != '' : sample += '_'
	x   = array('d', [0.5, 1.5, 2.5])
	exl = array('d', [0, 0, 0])
	exh = array('d', [0, 0, 0])
	y   = []
	eyl = []
	eyh = []
	for pdfset in pdfsets :
		y.append(pdfset[2])
		eyh.append(pdfset[3])
		eyl.append(pdfset[4])
	y   = array('d', y  )
	eyl = array('d', eyl)
	eyh = array('d', eyh)
	ttvStyle.ttvStyle()
	canvas = ROOT.TCanvas('canvas', 'canvas')
	graph = ROOT.TGraphAsymmErrors(3, x, y, exl, exh, eyl, eyh)
	histo = ROOT.TH1D('acceptances', 'acceptances', 3, 0., 3.)
	graph.Draw()
	graph.GetXaxis().SetRangeUser(0., 3.)
	graph.Draw()
	histo.Draw('goff')
	for i, pdfset in enumerate(pdfsets) :
		histo.SetBinContent(i+1, pdfset[2])
		histo.GetXaxis().SetBinLabel(i+1, pdfset[0])
	delta = graph.GetYaxis().GetXmax()-graph.GetYaxis().GetXmin()
	histo.SetMaximum(max(graph.GetYaxis().GetXmax() + 0.2*delta, central + 0.2*delta))
	histo.SetMinimum(min(graph.GetYaxis().GetXmin(), central - 0.2*delta))
	histo.GetYaxis().SetTitle('Acceptance')
	histo.GetXaxis().SetLabelFont(42)
	histo.GetXaxis().SetLabelSize(0.06)
	graph.SetMarkerStyle(20)
	graph.Draw('sameP')
	ROOT.TGaxis.SetMaxDigits(3)

	line = ROOT.TLine(0, central, 3, central)
	line.SetLineColor(ROOT.kRed)
	line.SetLineStyle(7)
	line.SetLineWidth(2)
	line.Draw()

	latex = ROOT.TLatex()
	latex.SetNDC()
	latex.SetTextFont(62)
	latex.SetTextSize(0.04)
	latex.SetTextAlign(13)
	x_pos = ROOT.gStyle.GetPadLeftMargin() + 0.03
	y_pos = 1. - ROOT.gStyle.GetPadTopMargin() - 0.03
	latex.DrawLatex(x_pos, y_pos, 'CMS Simulation')

	leg_entries = []
	leg_entries.append([line, 'CTEQ', 'l'])
	leg = draw_legend(leg_entries)
	leg.Draw()

	canvas.Update()
	canvas.Print('%spdf_acceptances.pdf' % sample)
#	raw_input('ok? ')


def draw_legend(entries) :
	leg = ROOT.TLegend()
	for entry in entries :
		leg.AddEntry(entry[0], ' ' + entry[1], entry[2])

	# set position
	width = 0.17
	x = 0.63
	y = 0.91
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


if __name__ == '__main__' :

	# ttW
	pdfsets = []
	pdfsets.append(('MSTW2008', '\\pdfMSTW{} ', 0.00428543, 0.00004291, 0.00008705))
	pdfsets.append(('NNPDF2.0', '\\pdfNNPDF{}', 0.00426639, 0.00002546, 0.00002425))
	pdfsets.append(('CT10'    , '\\pdfCT{}   ', 0.00426224, 0.00002749, 0.00003023))
	central = 0.00421983005211805
	make_table(central, pdfsets, 'ttw')
	make_plot (central, pdfsets, 'ttw')

	# WZ
	pdfsets = []
	pdfsets.append(('MSTW2008', '\\pdfMSTW{} ', 0.00008988, 0.00000184, 0.00000205))
	pdfsets.append(('NNPDF2.0', '\\pdfNNPDF{}', 0.00009225, 0.00000227, 0.00000222))
	pdfsets.append(('CT10'    , '\\pdfCT{}   ', 0.00009001, 0.00000319, 0.00000255))
	central = 0.0000857803883223637
	make_table(central, pdfsets, 'wz')
	make_plot (central, pdfsets, 'wz')
