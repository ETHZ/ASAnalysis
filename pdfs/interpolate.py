#include 'TH2.h'
#include <string>
#include <iostream>

## interpolate() implements a two step smearing. Its arguments are the histogram that should be interpolated
## and the direction in which to interpolate first (normally the direction in which the histogram changes 
## most slowly)
## Allowed types are SW (equivalently NE), SE (equivalently NW), NS, EW
## The second interpolation uses a 'Swiss Cross' average (non-zero N, S, E, W neighbors)
## TH2F* interpolate(const TH2F* hist, std::string firstInterpolationDirection)
## increases binning by factor of two and interpolates in specified direction
## TH2F* rebin(TH2F* hist, std::string firstInterpolationDirection)

def interpolate( hist,  firstInterpolationDirection):
	## TH2F *histCopy = (TH2F*)hist.Clone()
	histCopy = hist.Clone()
	
	## xStepPlus, xStepMinus, yStepPlus, yStepMinus
	if(firstInterpolationDirection=='SW' or firstInterpolationDirection=='NE' or firstInterpolationDirection=='Santa Fe'):
		xStepPlus = 1
		xStepMinus = -1
		yStepPlus = 1
		yStepMinus = -1
	
	elif(firstInterpolationDirection=='NW' or firstInterpolationDirection=='SE'):
		xStepPlus = -1
		xStepMinus = 1
		yStepPlus = 1
		yStepMinus = -1
	
	elif(firstInterpolationDirection=='N' or firstInterpolationDirection=='S' or firstInterpolationDirection=='NS' or firstInterpolationDirection=='SN'):
		xStepPlus = 0
		xStepMinus = 0
		yStepPlus = 1
		yStepMinus = -1
	
	elif(firstInterpolationDirection=='E' or firstInterpolationDirection=='W' or firstInterpolationDirection=='EW' or firstInterpolationDirection=='WE'):
		xStepPlus = 1
		xStepMinus = -1
		yStepPlus = 0
		yStepMinus = 0
	
	else:
		## to avoid uninitialized variable warnings
		xStepPlus=0
		xStepMinus=yStepPlus=yStepMinus=0
		print firstInterpolationDirection, 'is not an allowed smearing first interpolation direction.\n Allowed first interpolation directions are SW (equivalently NE), SE (equivalently NW), NS, EW'
	 	return 0
	
	
	## make temporary histograms to store the results of both steps
	hist_step1 = histCopy.Clone()
	hist_step1.Reset()
	hist_step2 = histCopy.Clone()
	hist_step2.Reset()
	
	nBinsX = histCopy.GetNbinsX()
	nBinsY = histCopy.GetNbinsY()
	
	foo = getHistMaxMinBins(histCopy)
	xMin, xMax, yMin, yMax = foo[0], foo[1], foo[2], foo[3]
	 
	for i in range(1, nBinsX+1):
		for j in range(1, nBinsY+1):
			## do not extrapolate outside the scan
			if(i<xMin or i>xMax or j<yMin or j>yMax or alongDiagonal(histCopy, i,j)): continue 
			binContent = histCopy.GetBinContent(i, j)
			binContentPlusStep = histCopy.GetBinContent(i+xStepPlus, j+yStepPlus)
			binContentMinusStep = histCopy.GetBinContent(i+xStepMinus, j+yStepMinus)
			nFilled = 0
			if(binContentPlusStep>0): nFilled+=1
			if(binContentMinusStep>0): nFilled+=1
			## if we are at an empty bin and there are neighbors
			## in specified direction with non-zero entries
			if(binContent==0 and nFilled>0):
				## average over non-zero entries
				binContent = (binContentPlusStep+binContentMinusStep)/nFilled
				hist_step1.SetBinContent(i,j,binContent)
	
	## add result of interpolation
	histCopy.Add(hist_step1)
	
	for i in range(1, nBinsX+1):
		for j in range(1, nBinsY+1):
			if(i<xMin or i>xMax or j<yMin or j>yMax or alongDiagonal(histCopy, i,j)): continue 
			binContent = histCopy.GetBinContent(i, j)
			## get entries for 'Swiss Cross' average
			binContentUp = histCopy.GetBinContent(i, j+1)
			binContentDown = histCopy.GetBinContent(i, j-1)
			binContentLeft = histCopy.GetBinContent(i-1, j)
			binContentRight = histCopy.GetBinContent(i+1, j)
			nFilled=0
			if(binContentUp>0): nFilled+=1
			if(binContentDown>0): nFilled+=1
			if(binContentRight>0): nFilled+=1
			if(binContentLeft>0): nFilled+=1
			if(binContent==0 and nFilled>0) :
				## only average over non-zero entries
				binContent = (binContentUp+binContentDown+binContentRight+binContentLeft)/nFilled
				hist_step2.SetBinContent(i,j,binContent)
		
	## add 'Swiss Cross' average
	histCopy.Add(hist_step2)
	
	return histCopy



## find absolute boundaries of the scan
## in most inefficient way possible (inefficient => simpler => hopefully fewer typos)
def getHistMaxMinBins(h):
	xMin=h.GetNbinsX() ## maximum possible minimum -- large dummy value
	yMin=h.GetNbinsY() ## large dummy value
	xMax=yMax=0
	
	for iX in range(1,h.GetNbinsX()+1):
		for iY in range(1,h.GetNbinsY()+1):
			if(h.GetBinContent(iX, iY)>1e-10):
				if(iX<xMin): xMin=iX
				if(iY<yMin): yMin=iY
				if(iX>xMax): xMax=iX
				if(iY>yMax): yMax=iY
	return xMin, xMax, yMin, yMax

## Omit bins with even bin index.
## Use for test in which known values are omitted
## to determine bias from interpolation
def omitEven(h):
	for i in range(0,h.GetNbinsX()+1):
		for j in range(0, jh.GetNbinsY()+1):
			if(i%2==0 or j%2==0): h.SetBinContent(i, j, 0)

## Omit bins with odd bin index.
## Use for test in which known values are omitted
## to determine bias from interpolation
def omitOdd(h):
	for i in range(0,h.GetNbinsX()+1):
		for  j in range(0, jh.GetNbinsY()+1):
			if((i+1)%2==0 or (j+1)%2==0): h.SetBinContent(i, j, 0)

## Tests if a bin is along a diagonal of the scan.
## Don't want to extrapolate past edge of scan in diagonal direction
def alongDiagonal(h, iX, iY):
	## calculate three most 'northwestern' neigbors
	sumNW = h.GetBinContent(iX, iY+1)+h.GetBinContent(iX-1, iY-1)+h.GetBinContent(iX-1, iY)
	## calculate three most 'southeastern' neigbors
	sumSE = h.GetBinContent(iX, iY-1)+h.GetBinContent(iX+1, iY-1)+h.GetBinContent(iX+1, iY)
	## etc.
	sumSW = h.GetBinContent(iX, iY-1)+h.GetBinContent(iX-1, iY-1)+h.GetBinContent(iX-1, iY)
	sumNE = h.GetBinContent(iX, iY+1)+h.GetBinContent(iX+1, iY+1)+h.GetBinContent(iX+1, iY)
	
	if((sumNW==0 and sumSE!=0) or (sumNW!=0 and sumSE==0) or (sumSW==0 and sumNE!=0) or (sumSW!=0 and sumNE==0)): return True
	else: return False

def rebin(hist, firstInterpolationDirection):
	histName = hist.GetName()
	histName+='_rebin'
	
	## bin widths are needed so as to not shift histogram by half a bin with each rebinning
	## assume constant binning
	binWidthX = hist.GetXaxis().GetBinWidth(1)
	binWidthY = hist.GetYaxis().GetBinWidth(1)
	
	histRebinned = ROOT.TH2F(histName, histName,
	  		  2*hist.GetNbinsX(),
	  		  hist.GetXaxis().GetXmin()+binWidthX/4,
	  		  hist.GetXaxis().GetXmax()+binWidthX/4,
	  		  2*hist.GetNbinsY(),
	  		  hist.GetYaxis().GetXmin()+binWidthY/4,
	  		  hist.GetYaxis().GetXmax()+binWidthY/4)
	
	## copy results from previous histogram
	for iX in range(1,h.GetNbinsX()+1):
		for iY in range(1,h.GetNbinsY()+1):
			binContent = hist.GetBinContent(iX, iY)
			histRebinned.SetBinContent(2*iX-1, 2*iY-1, binContent)
	histRebinned.SetMaximum(hist.GetMaximum())
	histRebinned.SetMinimum(hist.GetMinimum())
	
	## use interpolation to re-fill histogram
	histRebinnedInterpolated = interpolate(histRebinned, firstInterpolationDirection)
	
	return histRebinnedInterpolated
