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
pdfset = args[3]

global xvar, yvar

print 'running on model', model
currentDir = '/shome/lbaeni/workspace/ttW/CMSSW_5_3_7_patch5/src/ASAnalysis/pdfs/'


#if model == 'TTWJets':
if model == 'WZTo3LNu':
	## define input files
	#subdir1 = currentDir+'Nov06-MSTW2008nlo68cl-NNPDF20_100_TTWJets_madgraph/'
	#subdir2 = currentDir+'Nov06-MSTW2008nlo68cl_+68cl-MSTW2008nlo68cl_+68clhalf_TTWJets_madgraph/'
	#subdir3 = currentDir+'Nov06-MSTW2008nlo68cl_-68cl-MSTW2008nlo68cl_-68clhalf_TTWJets_madgraph/'
	#subdir4 = currentDir+'Nov06-NNPDF20_as_0116_100-NNPDF20_as_0117_100_TTWJets_madgraph/'
	#subdir5 = currentDir+'Nov11-NNPDF20_as_0118_100-NNPDF20_as_0120_100_TTWJets_madgraph/'
	#subdir6 = currentDir+'Nov06-NNPDF20_as_0121_100-NNPDF20_as_0122_100_TTWJets_madgraph/'
	#subdir7 = currentDir+'Nov11-cteq6mE-CT10-CT10as_TTWJets_madgraph/'
	subdir1 = currentDir+'Mar10-MSTW2008nlo68cl-NNPDF20_100_WZTo3LNu/'
	subdir2 = currentDir+'Mar10-MSTW2008nlo68cl_+68cl-MSTW2008nlo68cl_+68clhalf_WZTo3LNu/'
	subdir3 = currentDir+'Mar10-MSTW2008nlo68cl_-68cl-MSTW2008nlo68cl_-68clhalf_WZTo3LNu/'
	subdir4 = currentDir+'Mar10-NNPDF20_as_0116_100-NNPDF20_as_0117_100_WZTo3LNu/'
	subdir5 = currentDir+'Mar10-NNPDF20_as_0118_100-NNPDF20_as_0120_100_WZTo3LNu/'
	subdir6 = currentDir+'Mar10-NNPDF20_as_0121_100-NNPDF20_as_0122_100_WZTo3LNu/'
	subdir7 = currentDir+'Mar10-cteq6mE-CT10-CT10as_WZTo3LNu/'

	## set number of variations
	npdf = {}
	for subdir in [subdir1, subdir2, subdir3, subdir4, subdir5, subdir6, subdir7]:
		npdf[subdir] = {}
	npdf[subdir1]['pdf2'] = 40
	npdf[subdir1]['pdf3'] = 100
	npdf[subdir2]['pdf2'] = 40
	npdf[subdir2]['pdf3'] = 40
	npdf[subdir3]['pdf2'] = 40
	npdf[subdir3]['pdf3'] = 40
	npdf[subdir4]['pdf2'] = 5
	npdf[subdir4]['pdf3'] = 27
	npdf[subdir5]['pdf2'] = 72
	npdf[subdir5]['pdf3'] = 72
	npdf[subdir6]['pdf2'] = 27
	npdf[subdir6]['pdf3'] = 5
	npdf[subdir7]['pdf2'] = 52
	npdf[subdir7]['pdf3'] = 10

	if pdfset == 'CT10':
		subdirs = [subdir7]
		pdfsets = {}
		pdfsets[subdir7] = ['pdf2', 'pdf3']
		centralSubdir = subdir7
		centralPDF = 'pdf2'
		DeltaA = []
	elif pdfset == 'MSTW2008':
		subdirs = [subdir1, subdir2, subdir3]
		pdfsets = {}
		pdfsets[subdir1] = ['pdf2']
		pdfsets[subdir2] = ['pdf2', 'pdf3']
		pdfsets[subdir3] = ['pdf2', 'pdf3']
		centralSubdir = subdir1
		centralPDF = 'pdf2'
		Aalpha = []
	elif pdfset == 'NNPDF20':
		subdirs = [subdir1, subdir4, subdir5, subdir6]
		pdfsets = {}
		pdfsets[subdir1] = ['pdf3']
		pdfsets[subdir4] = ['pdf2', 'pdf3']
		pdfsets[subdir5] = ['pdf2', 'pdf3']
		pdfsets[subdir6] = ['pdf2', 'pdf3']
		centralSubdir = subdir1
		centralPDF = 'pdf3'
		errsNNPDF = [0., 0., 0., 0.]
	else:
		print 'SOMETHING WENT TERRIBLY WRONG!! CHECK THE PDF SET'
		sys.exit(0)
	xmin = 50
	xmax = 50
	ymin = 50
	ymax = 50
	step = 1
	xvar = 'foobar'
	yvar = 'foobar'
	zmax = 15.
#	modelRoot = 'TTWJets_madgraph'
	modelRoot = 'WZTo3LNu'
	relevantRegion = 'TTbarWSel'
	#relevantRegion = 'TTbarWPresel'
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
allHistos['pdf3'] = {}

for region in regions:
	if not region == relevantRegion: continue

	centralSubdir
	centralPDF
	centralFile = ROOT.TFile(centralSubdir+model+'_'+region+'.root', 'READ', 'centralFile')
	histo = centralFile.Get(modelRoot+'_nPass_'+centralPDF+'_0').Clone()
	histo.SetName(modelRoot+'_eff_'+centralPDF+'_0')
	histo.Divide(centralFile.Get(modelRoot+'_nTot_'+centralPDF+'_0'))
	bin = histo.FindBin(50.,50.)
	centralEff = histo.GetBinContent(bin)
	print 'A0: %.8f' %(centralEff)

	for subdir in subdirs:
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
		#resultFile2 = ROOT.TFile(subdir+model+'_TTbarWPresel.root', 'READ', 'resultFile2')
		resultFile.cd()

		##for pset in ['cteq', 'mstw', 'pdf1', 'cteq66']:
		##for pset in ['pdf1', 'pdf2']:
		for pset in pdfsets[subdir]:
			print pset
			print subdir
			allHistos[pset][region] = {}
			allHistos[pset][region] = {}
			for i in range(npdf[subdir][pset]):
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
					for i in range(npdf[subdir][pset]):
						binEff = allHistos[pset][region]['eff'+str(i)].GetBinContent(bin)
						binEffs.append(binEff)
						print '%2.0f binEff: %.8f' %(i, binEff)
					if pdfset == 'CT10':
						errs = helper.masterFormulaCT10(binEffs, centralEff)
						DeltaA.append([errs[0], errs[1]])
					elif pdfset == 'MSTW2008':
						print 'MSTW2008'
						errs = helper.masterFormula(binEffs)
						Aalpha.append([binEffs[0]+errs[0],binEffs[0]-errs[1]])
						print 'A + DeltaA+ = %.8f + %.8f = %.8f' %(binEffs[0], errs[0], binEffs[0]+errs[0])
						print 'A - DeltaA- = %.8f - %.8f = %.8f' %(binEffs[0], errs[0], binEffs[0]-errs[1])
					elif pdfset == 'NNPDF20':
						errs = helper.masterFormulaNNPDF(binEffs, centralEff)
						errsNNPDF[0] += errs[0]
						errsNNPDF[1] += errs[1]
						errsNNPDF[2] += errs[2]
						errsNNPDF[3] += errs[3]
						#print errs[0], errs[1], errs[2], errs[3]
						if (errs[0] != 1 and errs[1] != 1): print 'NNPDF2.0 errs: %.8f    %.8f' %(sqrt(errs[2]/(errs[0]-1.))/centralEff*100., sqrt(errs[3]/(errs[1]-1.))/centralEff*100.) 
					else:
						errs = helper.masterFormula(binEffs)
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
	if pdfset == 'CT10':
		print 'CT10 pdf:         DeltaA+ = %.8f   DeltaA- = %.8f' %(DeltaA[0][0], DeltaA[0][1])
		print 'CT10 alpha_s:     DeltaA+ = %.8f   DeltaA- = %.8f' %(DeltaA[1][0], DeltaA[1][1])
		print 'CT10 pdf+alpha_s: DeltaA+ = %.8f   DeltaA- = %.8f' %(sqrt(DeltaA[0][0]*DeltaA[0][0] + DeltaA[1][0]*DeltaA[1][0]), sqrt(DeltaA[0][1]*DeltaA[0][1] + DeltaA[1][1]*DeltaA[1][1]))
	if pdfset == 'NNPDF20':
		print errsNNPDF
		Ap = sqrt(errsNNPDF[2]/(errsNNPDF[0]-1.))
		Am = sqrt(errsNNPDF[3]/(errsNNPDF[1]-1.))
		print 'NNPDF20:   DeltaA+: %.8f   DeltaA-: %.8f' %(Ap, Am)
		print 'NNPDF20 errs: %.8f   %.8f' %(Ap/centralEff*100., Am/centralEff*100.)
	if pdfset == 'MSTW2008':
		AmaxP = 0.
		AmaxM = 0.
		AminM = 0.
		for alpha in range(len(Aalpha)):
			if AmaxP < Aalpha[alpha][0]: AmaxP = Aalpha[alpha][0]
			if AmaxM < Aalpha[alpha][1]: AmaxM = Aalpha[alpha][1]
			if (AminM > Aalpha[alpha][1] or AminM == 0) : AminM = Aalpha[alpha][1]
			print Aalpha[alpha]
		Ap = AmaxP - centralEff
		Am = centralEff - AminM
		print 'MSTW2008:   DeltaA+: %.8f   DeltaA-: %.8f' %(Ap, Am)
		print 'MSTW2008:   errs: %.8f   %.8f' %(Ap/centralEff*100., Am/centralEff*100.)
		print (AmaxP - centralEff, centralEff - AmaxM)

print 'DONE!!'
print '\n\n EVERYTHING is done... closing...'	
