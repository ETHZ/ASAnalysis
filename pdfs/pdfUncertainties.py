#! /bin/env python
import os, ROOT, helper, sys, copy
import interpolate
from math import *

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetNdivisions(310, 'x')
ROOT.gStyle.SetNdivisions(308, 'x')
ROOT.gStyle.SetNdivisions(308, 'y')

print 'This script produces the datacards for some regions and two models'

args = sys.argv
print args
model = args[1]

global xvar, yvar

print 'running on model', model
currentDir = '/shome/lbaeni/workspace/ttW/CMSSW_5_3_7_patch5/src/ASAnalysis/pdfs/'


if model == 'TTWJets':
	subdir  = currentDir+'Oct29-CT10-CT10as-MSTW2008nlo68cl_TTWJets_madgraph/'
	subdir2 = currentDir+'Oct29-MSTW2008nnlo68cl_asmz-68clhalf-NNPDF20_100_TTWJets_madgraph/'
	xmin = 50
	xmax = 50
	ymin = 50
	ymax = 50
	step = 1
	xvar = 'foobar'
	yvar = 'foobar'
	zmax = 15.
	modelRoot = 'TTWJets_madgraph'
	relevantRegion = 'TTbarWSel'
	modelPub = 'TTWJets'

else:
	print 'SOMETHING WENT TERRIBLY WRONG!! CHECK THE MODEL'
	sys.exit(0)


xrange = range(xmin, xmax+step, step)
yrange = range(ymin, ymax+step, step)

print xrange, yrange

## print xrange, yrange
## print len(xrange), len(yrange)

def setAxisTitle(histo):
	histo.GetXaxis().SetTitle(xvar)
	histo.GetYaxis().SetTitle(yvar)

whichRegions = args[2]
if whichRegions == 'ttw':
	regions = ['TTbarWSel', 'TTbarWPresel']
else:
	print 'SOMETHING WENT WRONG WITH SPECIFYING THE REGIONS'
	sys.exit(0)

allAcceptances    = {}

#ncteq, npdf1, nmstw, ncteq66 = 40, 52, 40, 44
npdf1, npdf2 = 40, 10

allHistos = {}
allHistos['pdf1'] = {}
allHistos['pdf2'] = {}

for region in regions:
	if not region == relevantRegion: continue
	infileName = subdir+model+'_'+region+'.root'
	if not os.path.isfile(infileName):
		print 'you\'re not reading a file, exiting...'
		sys.exit()
	canvasName  = 'foo'+region

	mystyle = helper.set1DStyle()
	mystyle.cd()

	print 'Running on region '+region
	
	print '... for the model'+model

	resultFile = ROOT.TFile(infileName, 'READ', 'resultFile')
	resultFile.cd()

	##for pset in ['cteq', 'mstw', 'pdf1', 'cteq66']:
	for pset in ['pdf1', 'pdf2']:
		print pset
		allHistos[pset][region] = {}
		allHistos[pset][region] = {}
		for i in range(npdf1):
			if pset == 'pdf2' and i > 9: continue
			allHistos[pset][region]['eff'+str(i)] = resultFile.Get(modelRoot+'_nPass_'+pset+'_'+str(i)).Clone()
			allHistos[pset][region]['eff'+str(i)].SetName(modelRoot+'_eff_'+pset+'_'+str(i))
			allHistos[pset][region]['eff'+str(i)].Divide(resultFile.Get(modelRoot+'_nTot_'+pset+'_' +str(i)))

		## pdfUncertainty  = ROOT.TH2D (region+'_pdfUncertainty_'+pset, region+'_pdfUncertainty_'+pset, len(xrange)-1, xmin, xmax, len(yrange)-1, ymin, ymax)
		pdfUncertainty  = ROOT.TH2D (region+'_pdfUncertainty_'+pset, region+'_pdfUncertainty_'+pset, 1, 50, 51, 1, 50, 51)

		for x in xrange:
			for y in yrange:
				print x, y
				bin = allHistos[pset][region]['eff1'].FindBin(x, y)
				binEffs = []
				for i in range(npdf1):
					if pset == 'pdf2' and i > 9: continue
					binEff = allHistos[pset][region]['eff'+str(i)].GetBinContent(bin)
					binEffs.append(binEff)
					print 'binEff: %.8f' %(binEff)
				errs = helper.masterFormulaCT10(binEffs)
				print 'errs: %.8f    %.8f' %(errs[0], errs[1]) 
				if binEff != 0.:
					##print '%.5f   %.5f   %.5f' %(binEff, errs[0]/binEff, errs[1]/binEff)
					pdfUncertainty.Fill(x, y, 100.*max(errs[0]/binEffs[0], errs[1]/binEffs[0]))
		canv = ROOT.TCanvas('foo', 'bar', 900, 675)
		helper.useNiceColorPalette()
		canv.cd()
		texts = helper.getModelText(model)
		pdfUncertainty.Draw('colz text')
		pdfUncertainty.GetZaxis().SetRangeUser(0.50, 1.50)
		pdfUncertainty.GetXaxis().SetTitle(xvar)
		pdfUncertainty.GetYaxis().SetTitle(yvar)
		helper.drawLatex(texts[0], 0.15, 0.82, 0.05)
		helper.drawLatex(texts[1], 0.15, 0.75, 0.05)
		if len(texts) > 2: helper.drawLatex(texts[2], 0.2, 0.68, 0.05)
		helper.drawLatex(pset+' PDF uncertainty (%) - '+region, 0.15, 0.94, 0.05)
		canv.SaveAs(model+'_pdfUncertianty_'+pset+'_'+region+'.pdf')
	

print 'DONE!!'
print '\n\n EVERYTHING is done... closing...'	
