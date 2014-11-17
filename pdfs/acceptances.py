#!/usr/bin/python
import ROOT
import os, sys
from array import array
import imp
ttvStyle = imp.load_source('ttvStyle', '../python/ttvStyle.py')
from math import log10
helper = imp.load_source('helper', '../python/helper.py')


def make_table(central, pdfsets, sample = '') :
	'''print latex table'''

	helper.mkdir('tables')
	exponent = int(log10(central))
	if exponent < 0 : exponent -= 1
	if sample != '' : sample += '_'
	table_name = 'tables/%spdf_acceptances.tex' % sample
	print '[status] writing %s' % table_name
	with open(table_name, 'w') as file :
		file.write('%!TEX root = ../../Dissertation.tex\n\n')
		file.write('\\sisetup{\n')
		file.write('\ttable-number-alignment = center,\n')
		file.write('\ttable-figures-integer  = 1,\n')
		file.write('\ttable-figures-decimal  = 2,\n')
		file.write('\ttable-figures-exponent = 1,\n')
		file.write('\ttable-sign-exponent    = true,\n')
		file.write('\tscientific-notation    = fixed,\n')
		file.write('\tfixed-exponent         = %d,\n' % exponent)
		file.write('\tround-mode             = places,\n')
		file.write('\tround-precision        = 2\n')
		file.write('}\n\n')
		file.write('\\begin{tabular}{lSSS}\n')
		file.write('\t\\toprule\n')
		file.write('\t{PDF set}   & {$A_0$} & {$\\Delta A^+$} & {$\\Delta A^-$} \\\\\n')
		file.write('\t\\midrule\n')
		file.write('\t\\pdfCTEQ{}  & %.2e & & \\\\\n' % central)
		file.write('\t\\midrule\n')
		for pdfset in pdfsets :
			file.write('\t%s & %.2e & %.2e & %.2e \\\\\n' % pdfset[1:])
		file.write('\t\\bottomrule\n')
		file.write('\\end{tabular}')


def make_plot(central, pdfsets, sample = '', TeX_switch = False) :
	helper.mkdir('figures')
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
	pl = ttvStyle.ttvStyle(lumi = -1, cms_label = 1, TeX_switch = TeX_switch)
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

	pl.draw_cmsLine()

	leg_entries = []
	leg_entries.append([line, 'CTEQ', 'l'])
	leg = pl.draw_legend(leg_entries)
	leg.Draw()

	canvas.Update()
	if TeX_switch :
		canvas.Print('figures/%spdf_acceptances.tex' % sample)
	else :
		canvas.Print('figures/%spdf_acceptances.pdf' % sample)
#	raw_input('ok? ')


if __name__ == '__main__' :
	args = sys.argv
	TeX_switch = False

	if ('--tex' in args) :
		TeX_switch = True

	# ttW
	pdfsets = []
	pdfsets.append(('CT10'    , '\\pdfCT{}   ', 0.00426224, 0.00002749, 0.00003023))
	pdfsets.append(('MSTW2008', '\\pdfMSTW{} ', 0.00428543, 0.00004291, 0.00008705))
	pdfsets.append(('NNPDF2.0', '\\pdfNNPDF{}', 0.00426639, 0.00002546, 0.00002425))
	central = 0.00421983005211805
	make_table(central, pdfsets, 'ttw')
	make_plot (central, pdfsets, 'ttw', TeX_switch)

	# WZ
	pdfsets = []
	pdfsets.append(('CT10'    , '\\pdfCT{}   ', 0.00009001, 0.00000319, 0.00000255))
	pdfsets.append(('MSTW2008', '\\pdfMSTW{} ', 0.00008988, 0.00000184, 0.00000205))
	pdfsets.append(('NNPDF2.0', '\\pdfNNPDF{}', 0.00009225, 0.00000227, 0.00000222))
	central = 0.0000857803883223637
	make_table(central, pdfsets, 'wz')
	make_plot (central, pdfsets, 'wz', TeX_switch)
