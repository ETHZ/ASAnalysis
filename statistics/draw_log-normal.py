#! /usr/bin/python

import ROOT
import imp
helper = imp.load_source('helper', '../python/helper.py')
ttvStyle = imp.load_source('ttvStyle', '../python/ttvStyle.py')


def get_logNormal(kappa, xmin = 0., xmax = 3.) :
	func = ROOT.TF1('func', 'ROOT::Math::lognormal_pdf(x, 0., TMath::Log(%f))' % kappa, xmin, xmax)
	return func

pl = ttvStyle.ttvStyle(lumi = -1, cms_label = 1, TeX_switch = True)
canvas = pl.get_canvas()
legend = []
ymax = 5.
canvas.DrawFrame(0., 0., 3., ymax, ';\\varepsilon = \\frac{\\theta}{\\tilde{\\theta}};\\text{Probability density } \\frac{\\mathrm{d}p}{\\mathrm{d}\\varepsilon}')
colors = []
colors.append(ROOT.kBlack)
colors.append(ROOT.kRed)
colors.append(ROOT.kGreen+4)
colors.append(ROOT.kBlue)

for i, kappa in enumerate([1.10, 1.20, 1.33, 1.50]) :
	func = get_logNormal(kappa)
	func.SetMaximum(ymax)
	func.SetLineColor(colors[i])
	func.SetLineStyle(i+1)
	print func.GetLineWidth()
	func.DrawCopy('C')
	legend.append([func, '\kappa = %4.2f' % kappa, 'l'])

leg = pl.draw_legend(legend)
leg.Draw()

print canvas.ls()

canvas.Update()
pl.save_canvas(canvas, './', 'log-normal')
#raw_input('ok?')
