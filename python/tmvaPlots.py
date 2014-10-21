#! /usr/bin/python
import ROOT
import os, sys
import helper
import ttvStyle


def plot_correlation(path, sb_switch = 'S') :
	if not os.path.exists(path) :
		print '[ERROR] %s does not exist!' % path
		sys.exit(1)

	pl = ttvStyle.ttvStyle(lumi = -1, cms_label = 1, TeX_switch = False)

	name = 'CorrelationMatrix%s' % sb_switch
	file = ROOT.TFile.Open(path, 'READ')
	histo = file.Get(name)

	canvas = pl.get_canvas()
	histo.SetMarkerColor(1)
#	histo.SetMarkerSize(1.4)
	histo.Scale(1./100.)
	histo.GetXaxis().LabelsOption('h')
	off = histo.GetYaxis().GetLabelOffset()
	histo.GetXaxis().SetLabelOffset(off)
	for x in range(1, histo.GetNbinsX()) :
		for y in range(x+1, histo.GetNbinsY()+1) :
			histo.SetBinContent(x, y, 0.)
#	histo.GetYaxis().SetRangeUser(0,6.2)
#	pl.ttvStyle.cd()
#	histo.GetXaxis().SetBinLabel(1, '\\mathrm{p_{\\perp}} (\\rm{\\ell})')
#	histo.GetXaxis().SetBinLabel(2, 'p_{T}(l_{2})')
#	histo.GetXaxis().SetBinLabel(4, 'E')
	histo.Draw('text')
	pl.draw_cmsLine()
	canvas.Update()
	canvas.Print('%s.pdf' % name)
#	raw_input('ok?')


if __name__ == '__main__' :
	args = sys.argv
	path = 'TMVA.root'

	if ('-f' in args) and (args[args.index('-f')+1] != '') :
		path = str(args[args.index('-f')+1])

	plot_correlation(path, 'S')
	plot_correlation(path, 'B')
